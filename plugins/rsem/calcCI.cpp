#include<ctime>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<algorithm>
#include<vector>
#include<pthread.h>

#include "utils.h"
#include "my_assert.h"
#include "sampling.h"

#include "Model.h"
#include "SingleModel.h"
#include "SingleQModel.h"
#include "PairedEndModel.h"
#include "PairedEndQModel.h"

#include "Refs.h"
#include "GroupInfo.h"

#include "Buffer.h"

using namespace std;

struct Params {
	int no;
	FILE *fi;
	engine_type *engine;
	const double *mw;
};

struct CIType {
	float lb, ub; // the interval is [lb, ub]

	CIType() { lb = ub = 0.0; }
};

struct CIParams {
	int no;
	int start_gene_id, end_gene_id;
};

int model_type;

int nMB;
double confidence;
int nCV, nSpC, nSamples; // nCV: number of count vectors; nSpC: number of theta vectors sampled per count vector; nSamples: nCV * nSpC
int nThreads;

float *l_bars;

char cvsF[STRLEN], tmpF[STRLEN], command[STRLEN];

CIType *iso_tpm, *gene_tpm, *iso_fpkm, *gene_fpkm;

int M, m;
Refs refs;
GroupInfo gi;
char imdName[STRLEN], statName[STRLEN];
char modelF[STRLEN], groupF[STRLEN], refF[STRLEN];

vector<double> eel; //expected effective lengths

Buffer *buffer;

bool quiet;

Params *paramsArray;
pthread_t *threads;
pthread_attr_t attr;
int rc;

CIParams *ciParamsArray;

template<class ModelType>
void calcExpectedEffectiveLengths(ModelType& model) {
	int lb, ub, span;
	double *pdf = NULL, *cdf = NULL, *clen = NULL; // clen[i] = \sigma_{j=1}^{i}pdf[i]*(lb+i)
  
	model.getGLD().copyTo(pdf, cdf, lb, ub, span);
	clen = new double[span + 1];
	clen[0] = 0.0;
	for (int i = 1; i <= span; i++) {
		clen[i] = clen[i - 1] + pdf[i] * (lb + i);
	}

	eel.assign(M + 1, 0.0);
	for (int i = 1; i <= M; i++) {
		int totLen = refs.getRef(i).getTotLen();
		int fullLen = refs.getRef(i).getFullLen();
		int pos1 = max(min(totLen - fullLen + 1, ub) - lb, 0);
		int pos2 = max(min(totLen, ub) - lb, 0);

		if (pos2 == 0) { eel[i] = 0.0; continue; }
    
		eel[i] = fullLen * cdf[pos1] + ((cdf[pos2] - cdf[pos1]) * (totLen + 1) - (clen[pos2] - clen[pos1]));
		assert(eel[i] >= 0);
		if (eel[i] < MINEEL) { eel[i] = 0.0; }
	}
  
	delete[] pdf;
	delete[] cdf;
	delete[] clen;
}

void* sample_theta_from_c(void* arg) {
	int *cvec;
	double *theta;
	float *tpm;
	gamma_dist **gammas;
	gamma_generator **rgs;

	Params *params = (Params*)arg;
	FILE *fi = params->fi;
	const double *mw = params->mw;

	cvec = new int[M + 1];
	theta = new double[M + 1];
	gammas = new gamma_dist*[M + 1];
	rgs = new gamma_generator*[M + 1];
	tpm = new float[M + 1];
	float l_bar; // the mean transcript length over the sample

	int cnt = 0;
	while (fscanf(fi, "%d", &cvec[0]) == 1) {
		for (int j = 1; j <= M; j++) assert(fscanf(fi, "%d", &cvec[j]) == 1);

		++cnt;

		for (int j = 0; j <= M; j++) {
			gammas[j] = new gamma_dist(cvec[j]);
			rgs[j] = new gamma_generator(*(params->engine), *gammas[j]);
		}

		for (int i = 0; i < nSpC; i++) {
			double sum = 0.0;
			for (int j = 0; j <= M; j++) {
				theta[j] = ((j == 0 || (eel[j] >= EPSILON && mw[j] >= EPSILON)) ? (*rgs[j])() / mw[j] : 0.0);
				sum += theta[j];
			}
			assert(sum >= EPSILON);
			for (int j = 0; j <= M; j++) theta[j] /= sum;

			sum = 0.0;
			tpm[0] = 0.0;
			for (int j = 1; j <= M; j++)
				if (eel[j] >= EPSILON) {
					tpm[j] = theta[j] / eel[j];
					sum += tpm[j];
				}
				else assert(theta[j] < EPSILON);
			assert(sum >= EPSILON);
			l_bar = 0.0; // store mean effective length of the sample
			for (int j = 1; j <= M; j++) { tpm[j] /= sum; l_bar += tpm[j] * eel[j]; tpm[j] *= 1e6; }
			buffer->write(l_bar, tpm + 1); // ommit the first element in tpm
		}

		for (int j = 0; j <= M; j++) {
			delete gammas[j];
			delete rgs[j];
		}

		if (verbose && cnt % 100 == 0) { printf("Thread %d, %d count vectors are processed!\n", params->no, cnt); }
	}

	delete[] cvec;
	delete[] theta;
	delete[] gammas;
	delete[] rgs;
	delete[] tpm;

	return NULL;
}

template<class ModelType>
void sample_theta_vectors_from_count_vectors() {
	ModelType model;
	model.read(modelF);
	calcExpectedEffectiveLengths<ModelType>(model);

	int num_threads = min(nThreads, nCV);

	buffer = new Buffer(nMB, nSamples, M, l_bars, tmpF);

	paramsArray = new Params[num_threads];
	threads = new pthread_t[num_threads];

	char inpF[STRLEN];
	for (int i = 0; i < num_threads; i++) {
		paramsArray[i].no = i;
		sprintf(inpF, "%s%d", cvsF, i);
		paramsArray[i].fi = fopen(inpF, "r");
		paramsArray[i].engine = engineFactory::new_engine();
		paramsArray[i].mw = model.getMW();
	}

	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for (int i = 0; i < num_threads; i++) {
		rc = pthread_create(&threads[i], &attr, &sample_theta_from_c, (void*)(&paramsArray[i]));
		pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) in sample_theta_vectors_from_count_vectors!");
	}
	for (int i = 0; i < num_threads; i++) {
		rc = pthread_join(threads[i], NULL);
		pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) in sample_theta_vectors_from_count_vectors!");
	}

	/* destroy attribute */
	pthread_attr_destroy(&attr);
	delete[] threads;

	for (int i = 0; i < num_threads; i++) {
		fclose(paramsArray[i].fi);
		delete paramsArray[i].engine;
	}
	delete[] paramsArray;

	delete buffer; // Must delete here, force the content left in the buffer be written into the disk

	if (verbose) { printf("Sampling is finished!\n"); }
}

void calcCI(int nSamples, float *samples, float &lb, float &ub) {
	int p, q; // p pointer for lb, q pointer for ub;
	int newp, newq;
	int threshold = nSamples - (int(confidence * nSamples - 1e-8) + 1);
	int nOutside = 0;

	sort(samples, samples + nSamples);

	p = 0; q = nSamples - 1;
	newq = nSamples - 1;
	do {
		q = newq;
		while (newq > 0 && samples[newq - 1] == samples[newq]) newq--;
		newq--;
	} while (newq >= 0 && nSamples - (newq + 1) <= threshold);

	nOutside = nSamples - (q + 1);

	lb = -1e30; ub = 1e30;
	do {
		if (samples[q] - samples[p] < ub - lb) {
			lb = samples[p];
			ub = samples[q];
		}

		newp = p;
		while (newp < nSamples - 1 && samples[newp] == samples[newp + 1]) newp++;
		newp++;
		if (newp <= threshold) {
			nOutside += newp - p;
			p = newp;
			while (nOutside > threshold && q < nSamples - 1) {
				newq = q + 1;
				while (newq < nSamples - 1 && samples[newq] == samples[newq + 1]) newq++;
				nOutside -= newq - q;
				q = newq;
			}
			assert(nOutside <= threshold);
		}
		else p = newp;
	} while (p <= threshold);
}

void* calcCI_batch(void* arg) {
	float *itsamples, *gtsamples, *ifsamples, *gfsamples;
	ifstream fin;
	CIParams *ciParams = (CIParams*)arg;

	itsamples = new float[nSamples];
	gtsamples = new float[nSamples];
	ifsamples = new float[nSamples];
	gfsamples = new float[nSamples];

	fin.open(tmpF, ios::binary);
	// minus 1 here for that theta0 is not written!
	streampos pos = streampos(gi.spAt(ciParams->start_gene_id) - 1) * nSamples * FLOATSIZE;
	fin.seekg(pos, ios::beg);

	int cnt = 0;
	for (int i = ciParams->start_gene_id; i < ciParams->end_gene_id; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		memset(gtsamples, 0, FLOATSIZE * nSamples);
		memset(gfsamples, 0, FLOATSIZE * nSamples);
		for (int j = b; j < e; j++) {
			for (int k = 0; k < nSamples; k++) {
				fin.read((char*)(&itsamples[k]), FLOATSIZE);
				gtsamples[k] += itsamples[k];
				ifsamples[k] = 1e3 / l_bars[k] * itsamples[k];
				gfsamples[k] += ifsamples[k];
			}
			calcCI(nSamples, itsamples, iso_tpm[j].lb, iso_tpm[j].ub);
			calcCI(nSamples, ifsamples, iso_fpkm[j].lb, iso_fpkm[j].ub);
		}

		if (e - b > 1) {
			calcCI(nSamples, gtsamples, gene_tpm[i].lb, gene_tpm[i].ub);
			calcCI(nSamples, gfsamples, gene_fpkm[i].lb, gene_fpkm[i].ub);
		}
		else {
			gene_tpm[i].lb = iso_tpm[b].lb; gene_tpm[i].ub = iso_tpm[b].ub;
			gene_fpkm[i].lb = iso_fpkm[b].lb; gene_fpkm[i].ub = iso_fpkm[b].ub;
		}

		++cnt;
		if (verbose && cnt % 1000 == 0) { printf("In thread %d, %d genes are processed for CI calculation!\n", ciParams->no, cnt); }
	}

	fin.close();

	delete[] itsamples;
	delete[] gtsamples;

	return NULL;
}

void calculate_credibility_intervals(char* imdName) {
	FILE *fo;
	char outF[STRLEN];
	int num_threads = nThreads;

	iso_tpm = new CIType[M + 1];
	gene_tpm = new CIType[m];
	iso_fpkm = new CIType[M + 1];
	gene_fpkm = new CIType[m];

	assert(M > 0);
	int quotient = M / num_threads;
	if (quotient < 1) { num_threads = M; quotient = 1; }
	int cur_gene_id = 0;
	int num_isoforms = 0;

	// A just so so strategy for paralleling
	ciParamsArray = new CIParams[num_threads];
	for (int i = 0; i < num_threads; i++) {
		ciParamsArray[i].no = i;
		ciParamsArray[i].start_gene_id = cur_gene_id;
		num_isoforms = 0;

		while ((m - cur_gene_id > num_threads - i - 1) && (i == num_threads - 1 || num_isoforms < quotient)) {
			num_isoforms += gi.spAt(cur_gene_id + 1) - gi.spAt(cur_gene_id);
			++cur_gene_id;
		}

		ciParamsArray[i].end_gene_id = cur_gene_id;
	}

	threads = new pthread_t[num_threads];

	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	// paralleling
	for (int i = 0; i < num_threads; i++) {
		rc = pthread_create(&threads[i], &attr, &calcCI_batch, (void*)(&ciParamsArray[i]));
		pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) in calculate_credibility_intervals!");
	}
	for (int i = 0; i < num_threads; i++) {
		rc = pthread_join(threads[i], NULL);
		pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) in calculate_credibility_intervals!");
	}

	// releasing resources

	/* destroy attribute */
	pthread_attr_destroy(&attr);
	delete[] threads;

	delete[] ciParamsArray;

	//isoform level results
	sprintf(outF, "%s.iso_res", imdName);
	fo = fopen(outF, "a");
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.6g%c", iso_tpm[i].lb, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.6g%c", iso_tpm[i].ub, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.6g%c", iso_fpkm[i].lb, (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.6g%c", iso_fpkm[i].ub, (i < M ? '\t' : '\n'));
	fclose(fo);

	//gene level results
	sprintf(outF, "%s.gene_res", imdName);
	fo = fopen(outF, "a");
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.6g%c", gene_tpm[i].lb, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.6g%c", gene_tpm[i].ub, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.6g%c", gene_fpkm[i].lb, (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.6g%c", gene_fpkm[i].ub, (i < m - 1 ? '\t' : '\n'));
	fclose(fo);

	delete[] iso_tpm;
	delete[] gene_tpm;
	delete[] iso_fpkm;
	delete[] gene_fpkm;

	if (verbose) { printf("All credibility intervals are calculated!\n"); }
}

int main(int argc, char* argv[]) {
	if (argc < 8) {
		printf("Usage: rsem-calculate-credibility-intervals reference_name imdName statName confidence nCV nSpC nMB [-p #Threads] [-q]\n");
		exit(-1);
	}

	strcpy(imdName, argv[2]);
	strcpy(statName, argv[3]);

	confidence = atof(argv[4]);
	nCV = atoi(argv[5]);
	nSpC = atoi(argv[6]);
	nMB = atoi(argv[7]);

	nThreads = 1;
	quiet = false;
	for (int i = 8; i < argc; i++) {
		if (!strcmp(argv[i], "-p")) nThreads = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "-q")) quiet = true;
	}
	verbose = !quiet;

	sprintf(refF, "%s.seq", argv[1]);
	refs.loadRefs(refF, 1);
	M = refs.getM();
	sprintf(groupF, "%s.grp", argv[1]);
	gi.load(groupF);
	m = gi.getm();

	nSamples = nCV * nSpC;
	assert(nSamples > 0 && M > 0); // for Buffter.h: (bufsize_type)nSamples
	l_bars = new float[nSamples];

	sprintf(tmpF, "%s.tmp", imdName);
	sprintf(cvsF, "%s.countvectors", imdName);

	sprintf(modelF, "%s.model", statName);
	FILE *fi = fopen(modelF, "r");
	general_assert(fi != NULL, "Cannot open " + cstrtos(modelF) + "!");
	assert(fscanf(fi, "%d", &model_type) == 1);
	fclose(fi);

	// Phase I
	switch(model_type) {
	case 0 : sample_theta_vectors_from_count_vectors<SingleModel>(); break;
	case 1 : sample_theta_vectors_from_count_vectors<SingleQModel>(); break;
	case 2 : sample_theta_vectors_from_count_vectors<PairedEndModel>(); break;
	case 3 : sample_theta_vectors_from_count_vectors<PairedEndQModel>(); break;
	}

	// Phase II
	calculate_credibility_intervals(imdName);

	delete l_bars;
	/*
	sprintf(command, "rm -f %s", tmpF);
	int status = system(command);
	if (status != 0) {
		fprintf(stderr, "Cannot delete %s!\n", tmpF);
		exit(-1);
	}
	*/

	return 0;
}
