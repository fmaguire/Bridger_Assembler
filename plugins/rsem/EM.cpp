#include<ctime>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cassert>
#include<string>
#include<vector>
#include<algorithm>
#include<fstream>
#include<iostream>
#include<pthread.h>

#include "utils.h"
#include "my_assert.h"
#include "sampling.h"

#include "Read.h"
#include "SingleRead.h"
#include "SingleReadQ.h"
#include "PairedEndRead.h"
#include "PairedEndReadQ.h"

#include "SingleHit.h"
#include "PairedEndHit.h"

#include "Model.h"
#include "SingleModel.h"
#include "SingleQModel.h"
#include "PairedEndModel.h"
#include "PairedEndQModel.h"

#include "Transcript.h"
#include "Transcripts.h"

#include "Refs.h"
#include "GroupInfo.h"
#include "HitContainer.h"
#include "ReadIndex.h"
#include "ReadReader.h"

#include "ModelParams.h"

#include "HitWrapper.h"
#include "BamWriter.h"

using namespace std;

const double STOP_CRITERIA = 0.001;
const int MAX_ROUND = 10000;
const int MIN_ROUND = 20;

struct Params {
	void *model;
	void *reader, *hitv, *ncpv, *mhp, *countv;
};

int read_type;
int m, M; // m genes, M isoforms
READ_INT_TYPE N0, N1, N2, N_tot;
int nThreads;


bool genBamF; // If user wants to generate bam file, true; otherwise, false.
bool bamSampling; // true if sampling from read posterior distribution when bam file is generated
bool updateModel, calcExpectedWeights;
bool genGibbsOut; // generate file for Gibbs sampler

char refName[STRLEN], outName[STRLEN];
char imdName[STRLEN], statName[STRLEN];
char refF[STRLEN], groupF[STRLEN], cntF[STRLEN], tiF[STRLEN];
char mparamsF[STRLEN];
char modelF[STRLEN], thetaF[STRLEN];

char inpSamType;
char *pt_fn_list, *pt_chr_list;
char inpSamF[STRLEN], outBamF[STRLEN], fn_list[STRLEN], chr_list[STRLEN];

char out_for_gibbs_F[STRLEN];

vector<double> theta, eel; // eel : expected effective length

double *probv, **countvs;

Refs refs;
GroupInfo gi;
Transcripts transcripts;

ModelParams mparams;

template<class ReadType, class HitType, class ModelType>
void init(ReadReader<ReadType> **&readers, HitContainer<HitType> **&hitvs, double **&ncpvs, ModelType **&mhps) {
	READ_INT_TYPE nReads;
	HIT_INT_TYPE nHits;
	int rt; // read type

	READ_INT_TYPE nrLeft, curnr; // nrLeft : number of reads left, curnr: current number of reads
	HIT_INT_TYPE nhT; // nhT : hit threshold per thread
	char datF[STRLEN];

	int s;
	char readFs[2][STRLEN];
	ReadIndex *indices[2];
	ifstream fin;

	readers = new ReadReader<ReadType>*[nThreads];
	genReadFileNames(imdName, 1, read_type, s, readFs);
	for (int i = 0; i < s; i++) {
		indices[i] = new ReadIndex(readFs[i]);
	}
	for (int i = 0; i < nThreads; i++) {
		readers[i] = new ReadReader<ReadType>(s, readFs, refs.hasPolyA(), mparams.seedLen); // allow calculation of calc_lq() function
		readers[i]->setIndices(indices);
	}

	hitvs = new HitContainer<HitType>*[nThreads];
	for (int i = 0; i < nThreads; i++) {
		hitvs[i] = new HitContainer<HitType>();
	}

	sprintf(datF, "%s.dat", imdName);
	fin.open(datF);
	general_assert(fin.is_open(), "Cannot open " + cstrtos(datF) + "! It may not exist.");
	fin>>nReads>>nHits>>rt;
	general_assert(nReads == N1, "Number of alignable reads does not match!");
	general_assert(rt == read_type, "Data file (.dat) does not have the right read type!");


	//A just so so strategy for paralleling
	nhT = nHits / nThreads;
	nrLeft = N1;
	curnr = 0;

	ncpvs = new double*[nThreads];
	for (int i = 0; i < nThreads; i++) {
		HIT_INT_TYPE ntLeft = nThreads - i - 1; // # of threads left

		general_assert(readers[i]->locate(curnr), "Read indices files do not match!");

		while (nrLeft > ntLeft && (i == nThreads - 1 || hitvs[i]->getNHits() < nhT)) {
			general_assert(hitvs[i]->read(fin), "Cannot read alignments from .dat file!");

			--nrLeft;
			if (verbose && nrLeft % 1000000 == 0) { cout<< "DAT "<< nrLeft << " reads left"<< endl; }
		}
		ncpvs[i] = new double[hitvs[i]->getN()];
		memset(ncpvs[i], 0, sizeof(double) * hitvs[i]->getN());
		curnr += hitvs[i]->getN();

		if (verbose) { cout<<"Thread "<< i<< " : N = "<< hitvs[i]->getN()<< ", NHit = "<< hitvs[i]->getNHits()<< endl; }
	}

	fin.close();

	mhps = new ModelType*[nThreads];
	for (int i = 0; i < nThreads; i++) {
		mhps[i] = new ModelType(mparams, false); // just model helper
	}

	probv = new double[M + 1];
	countvs = new double*[nThreads];
	for (int i = 0; i < nThreads; i++) {
		countvs[i] = new double[M + 1];
	}


	if (verbose) { printf("EM_init finished!\n"); }
}

template<class ReadType, class HitType, class ModelType>
void* E_STEP(void* arg) {
	Params *params = (Params*)arg;
	ModelType *model = (ModelType*)(params->model);
	ReadReader<ReadType> *reader = (ReadReader<ReadType>*)(params->reader);
	HitContainer<HitType> *hitv = (HitContainer<HitType>*)(params->hitv);
	double *ncpv = (double*)(params->ncpv);
	ModelType *mhp = (ModelType*)(params->mhp);
	double *countv = (double*)(params->countv);

	bool needCalcConPrb = model->getNeedCalcConPrb();

	ReadType read;

	READ_INT_TYPE N = hitv->getN();
	double sum;
	vector<double> fracs; //to remove this, do calculation twice
	HIT_INT_TYPE fr, to, id;

	if (needCalcConPrb || updateModel) { reader->reset(); }
	if (updateModel) { mhp->init(); }

	memset(countv, 0, sizeof(double) * (M + 1));
	for (READ_INT_TYPE i = 0; i < N; i++) {
		if (needCalcConPrb || updateModel) {
			general_assert(reader->next(read), "Can not load a read!");
		}

		fr = hitv->getSAt(i);
		to = hitv->getSAt(i + 1);
		fracs.resize(to - fr + 1);

		sum = 0.0;

		if (needCalcConPrb) { ncpv[i] = model->getNoiseConPrb(read); }
		fracs[0] = probv[0] * ncpv[i];
		if (fracs[0] < EPSILON) fracs[0] = 0.0;
		sum += fracs[0];
		for (HIT_INT_TYPE j = fr; j < to; j++) {
			HitType &hit = hitv->getHitAt(j);
			if (needCalcConPrb) { hit.setConPrb(model->getConPrb(read, hit)); }
			id = j - fr + 1;
			fracs[id] = probv[hit.getSid()] * hit.getConPrb();
			if (fracs[id] < EPSILON) fracs[id] = 0.0;
			sum += fracs[id];
		}

		if (sum >= EPSILON) {
			fracs[0] /= sum;
			countv[0] += fracs[0];
			if (updateModel) { mhp->updateNoise(read, fracs[0]); }
			if (calcExpectedWeights) { ncpv[i] = fracs[0]; }
			for (HIT_INT_TYPE j = fr; j < to; j++) {
				HitType &hit = hitv->getHitAt(j);
				id = j - fr + 1;
				fracs[id] /= sum;
				countv[hit.getSid()] += fracs[id];
				if (updateModel) { mhp->update(read, hit, fracs[id]); }
				if (calcExpectedWeights) { hit.setConPrb(fracs[id]); }
			}			
		}
		else if (calcExpectedWeights) {
			ncpv[i] = 0.0;
			for (HIT_INT_TYPE j = fr; j < to; j++) {
				HitType &hit = hitv->getHitAt(j);
				hit.setConPrb(0.0);
			}
		}
	}

	return NULL;
}

template<class ReadType, class HitType, class ModelType>
void* calcConProbs(void* arg) {
	Params *params = (Params*)arg;
	ModelType *model = (ModelType*)(params->model);
	ReadReader<ReadType> *reader = (ReadReader<ReadType>*)(params->reader);
	HitContainer<HitType> *hitv = (HitContainer<HitType>*)(params->hitv);
	double *ncpv = (double*)(params->ncpv);

	ReadType read;
	READ_INT_TYPE N = hitv->getN();
	HIT_INT_TYPE fr, to;

	assert(model->getNeedCalcConPrb());
	reader->reset();

	for (READ_INT_TYPE i = 0; i < N; i++) {
		general_assert(reader->next(read), "Can not load a read!");

		fr = hitv->getSAt(i);
		to = hitv->getSAt(i + 1);

		ncpv[i] = model->getNoiseConPrb(read);
		for (HIT_INT_TYPE j = fr; j < to; j++) {
			HitType &hit = hitv->getHitAt(j);
			hit.setConPrb(model->getConPrb(read, hit));
		}
	}

	return NULL;
}

template<class ModelType>
void calcExpectedEffectiveLengths(ModelType& model) {
	int lb, ub, span;
	double *pdf = NULL, *cdf = NULL, *clen = NULL; // clen[i] = sigma_{j=1}^{i}pdf[i]*(lb+i)
  
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

void polishTheta(vector<double>& theta, const vector<double>& eel, const double* mw) {
	double sum = 0.0;

	/* The reason that for noise gene, mw value is 1 is :
	 * currently, all masked positions are for poly(A) sites, which in theory should be filtered out.
	 * So the theta0 does not containing reads from any masked position
	 */

	for (int i = 0; i <= M; i++) {
		// i == 0, mw[i] == 1
		if (i > 0 && (mw[i] < EPSILON || eel[i] < EPSILON)) {
			theta[i] = 0.0;
			continue;
		}
		theta[i] = theta[i] / mw[i];
		sum += theta[i];
	}
	// currently is OK, since no transcript should be masked totally, only the poly(A) tail related part will be masked
	general_assert(sum >= EPSILON, "No effective length is no less than" + ftos(MINEEL, 6) + " !");
	for (int i = 0; i <= M; i++) theta[i] /= sum;
}

void calcExpressionValues(const vector<double>& theta, const vector<double>& eel, vector<double>& tpm, vector<double>& fpkm) {
	double denom;
	vector<double> frac;

	//calculate fraction of count over all mappabile reads
	denom = 0.0;
	frac.assign(M + 1, 0.0);
	for (int i = 1; i <= M; i++) 
	  if (eel[i] >= EPSILON) {
	    frac[i] = theta[i];
	    denom += frac[i];
	  }
	general_assert(denom >= EPSILON, "No alignable reads?!");
	for (int i = 1; i <= M; i++) frac[i] /= denom;
  
	//calculate FPKM
	fpkm.assign(M + 1, 0.0);
	for (int i = 1; i <= M; i++)
		if (eel[i] >= EPSILON) fpkm[i] = frac[i] * 1e9 / eel[i];

	//calculate TPM
	tpm.assign(M + 1, 0.0);
	denom = 0.0;
	for (int i = 1; i <= M; i++) denom += fpkm[i];
	for (int i = 1; i <= M; i++) tpm[i] = fpkm[i] / denom * 1e6;  
}

template<class ModelType>
void writeResults(ModelType& model, double* counts) {
	char outF[STRLEN];
	FILE *fo;

	sprintf(modelF, "%s.model", statName);
	model.write(modelF);

	vector<int> tlens;
	vector<double> fpkm, tpm, isopct;
	vector<double> glens, gene_eels, gene_counts, gene_tpm, gene_fpkm;

	calcExpressionValues(theta, eel, tpm, fpkm);

	//calculate IsoPct, etc.
	isopct.assign(M + 1, 0.0);
	tlens.assign(M + 1, 0);

	glens.assign(m, 0.0); gene_eels.assign(m, 0.0);
	gene_counts.assign(m, 0.0); gene_tpm.assign(m, 0.0); gene_fpkm.assign(m, 0.0);

	for (int i = 0; i < m; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		for (int j = b; j < e; j++) {
			const Transcript& transcript = transcripts.getTranscriptAt(j);
			tlens[j] = transcript.getLength();

			gene_counts[i] += counts[j];
			gene_tpm[i] += tpm[j];
			gene_fpkm[i] += fpkm[j];
		}

		if (gene_tpm[i] < EPSILON) {
			double frac = 1.0 / (e - b);
			for (int j = b; j < e; j++) {
				glens[i] += tlens[j] * frac;
				gene_eels[i] += eel[j] * frac;
			}
		}
		else {
			for (int j = b; j < e; j++) {
				isopct[j] = tpm[j] / gene_tpm[i];
				glens[i] += tlens[j] * isopct[j];
				gene_eels[i] += eel[j] * isopct[j];
			}
		}
	}

	//isoform level results
	sprintf(outF, "%s.iso_res", imdName);
	fo = fopen(outF, "w");
	for (int i = 1; i <= M; i++) {
		const Transcript& transcript = transcripts.getTranscriptAt(i);
		fprintf(fo, "%s%c", transcript.getTranscriptID().c_str(), (i < M ? '\t' : '\n'));
	}
	for (int i = 1; i <= M; i++) {
		const Transcript& transcript = transcripts.getTranscriptAt(i);
		fprintf(fo, "%s%c", transcript.getGeneID().c_str(), (i < M ? '\t' : '\n'));
	}
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%d%c", tlens[i], (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", eel[i], (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", counts[i], (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", tpm[i], (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", fpkm[i], (i < M ? '\t' : '\n'));
	for (int i = 1; i <= M; i++)
		fprintf(fo, "%.2f%c", isopct[i] * 1e2, (i < M ? '\t' : '\n'));
	fclose(fo);

	//gene level results
	sprintf(outF, "%s.gene_res", imdName);
	fo = fopen(outF, "w");
	for (int i = 0; i < m; i++) {
		const Transcript& transcript = transcripts.getTranscriptAt(gi.spAt(i));
		fprintf(fo, "%s%c", transcript.getGeneID().c_str(), (i < m - 1 ? '\t' : '\n'));
	}
	for (int i = 0; i < m; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		for (int j = b; j < e; j++) {
			fprintf(fo, "%s%c", transcripts.getTranscriptAt(j).getTranscriptID().c_str(), (j < e - 1 ? ',' : (i < m - 1 ? '\t' :'\n')));
		}
	}
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", glens[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_eels[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_counts[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_tpm[i], (i < m - 1 ? '\t' : '\n'));
	for (int i = 0; i < m; i++)
		fprintf(fo, "%.2f%c", gene_fpkm[i], (i < m - 1 ? '\t' : '\n'));
	fclose(fo);

	if (verbose) { printf("Expression Results are written!\n"); }
}

template<class ReadType, class HitType, class ModelType>
void release(ReadReader<ReadType> **readers, HitContainer<HitType> **hitvs, double **ncpvs, ModelType **mhps) {
	delete[] probv;
	for (int i = 0; i < nThreads; i++) {
		delete[] countvs[i];
	}
	delete[] countvs;

	for (int i = 0; i < nThreads; i++) {
		delete readers[i];
		delete hitvs[i];
		delete[] ncpvs[i];
		delete mhps[i];
	}
	delete[] readers;
	delete[] hitvs;
	delete[] ncpvs;
	delete[] mhps;
}

inline bool doesUpdateModel(int ROUND) {
  //  return ROUND <= 20 || ROUND % 100 == 0;
  return ROUND <= 10;
}

//Including initialize, algorithm and results saving
template<class ReadType, class HitType, class ModelType>
void EM() {
	FILE *fo;

	int ROUND;
	double sum;

	double bChange = 0.0, change = 0.0; // bChange : biggest change
	int totNum = 0;

	ModelType model(mparams); //master model
	ReadReader<ReadType> **readers;
	HitContainer<HitType> **hitvs;
	double **ncpvs;
	ModelType **mhps; //model helpers

	Params fparams[nThreads];
	pthread_t threads[nThreads];
	pthread_attr_t attr;
	int rc;


	//initialize boolean variables
	updateModel = calcExpectedWeights = false;

	theta.clear();
	theta.resize(M + 1, 0.0);
	init<ReadType, HitType, ModelType>(readers, hitvs, ncpvs, mhps);

	//set initial parameters
	assert(N_tot > N2);
	theta[0] = max(N0 * 1.0 / (N_tot - N2), 1e-8);
	double val = (1.0 - theta[0]) / M;
	for (int i = 1; i <= M; i++) theta[i] = val;

	model.estimateFromReads(imdName);

	for (int i = 0; i < nThreads; i++) {
		fparams[i].model = (void*)(&model);

		fparams[i].reader = (void*)readers[i];
		fparams[i].hitv = (void*)hitvs[i];
		fparams[i].ncpv = (void*)ncpvs[i];
		fparams[i].mhp = (void*)mhps[i];
		fparams[i].countv = (void*)countvs[i];
	}

	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	ROUND = 0;
	do {
		++ROUND;

		updateModel = doesUpdateModel(ROUND);

		for (int i = 0; i <= M; i++) probv[i] = theta[i];

		//E step
		for (int i = 0; i < nThreads; i++) {
			rc = pthread_create(&threads[i], &attr, E_STEP<ReadType, HitType, ModelType>, (void*)(&fparams[i]));
			pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) at ROUND " + itos(ROUND) + "!");
		}

		for (int i = 0; i < nThreads; i++) {
			rc = pthread_join(threads[i], NULL);
			pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) at ROUND " + itos(ROUND) + "!");
		}

		model.setNeedCalcConPrb(false);

		for (int i = 1; i < nThreads; i++) {
			for (int j = 0; j <= M; j++) {
				countvs[0][j] += countvs[i][j];
			}
		}

		//add N0 noise reads
		countvs[0][0] += N0;

		//M step;
		sum = 0.0;
		for (int i = 0; i <= M; i++) sum += countvs[0][i];
		assert(sum >= EPSILON);
		for (int i = 0; i <= M; i++) theta[i] = countvs[0][i] / sum;

		if (updateModel) {
			model.init();
			for (int i = 0; i < nThreads; i++) { model.collect(*mhps[i]); }
			model.finish();
		}

		// Relative error
		bChange = 0.0; totNum = 0;
		for (int i = 0; i <= M; i++)
			if (probv[i] >= 1e-7) {
				change = fabs(theta[i] - probv[i]) / probv[i];
				if (change >= STOP_CRITERIA) ++totNum;
				if (bChange < change) bChange = change;
			}

		if (verbose) { cout<< "ROUND = "<< ROUND<< ", SUM = "<< setprecision(15)<< sum<< ", bChange = " << setprecision(6)<< bChange<< ", totNum = " << totNum<< endl; }
	} while (ROUND < MIN_ROUND || (totNum > 0 && ROUND < MAX_ROUND));
//	} while (ROUND < 1);

	if (totNum > 0) { cout<< "Warning: RSEM reaches "<< MAX_ROUND<< " iterations before meeting the convergence criteria."<< endl; }

	//generate output file used by Gibbs sampler
	if (genGibbsOut) {
		if (model.getNeedCalcConPrb()) {
			for (int i = 0; i < nThreads; i++) {
				rc = pthread_create(&threads[i], &attr, calcConProbs<ReadType, HitType, ModelType>, (void*)(&fparams[i]));
				pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) when generating files for Gibbs sampler!");
			}
			for (int i = 0; i < nThreads; i++) {
				rc = pthread_join(threads[i], NULL);
				pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) when generating files for Gibbs sampler!");
			}
		}
		model.setNeedCalcConPrb(false);

		sprintf(out_for_gibbs_F, "%s.ofg", imdName);
		ofstream fout(out_for_gibbs_F);
		fout<< M<< " "<< N0<< endl;
		for (int i = 0; i < nThreads; i++) {
			READ_INT_TYPE numN = hitvs[i]->getN();
			for (READ_INT_TYPE j = 0; j < numN; j++) {
				HIT_INT_TYPE fr = hitvs[i]->getSAt(j);
				HIT_INT_TYPE to = hitvs[i]->getSAt(j + 1);
				HIT_INT_TYPE totNum = 0;

				if (ncpvs[i][j] >= EPSILON) { ++totNum; fout<< "0 "<< setprecision(15)<< ncpvs[i][j]<< " "; }
				for (HIT_INT_TYPE k = fr; k < to; k++) {
					HitType &hit = hitvs[i]->getHitAt(k);
					if (hit.getConPrb() >= EPSILON) {
						++totNum;
						fout<< hit.getSid()<< " "<< setprecision(15)<< hit.getConPrb()<< " ";
					}
				}

				if (totNum > 0) { fout<< endl; }
			}
		}
		fout.close();
	}

	//calculate expected weights and counts using learned parameters
	//just use the raw theta learned from the data, do not correct for eel or mw
	updateModel = false; calcExpectedWeights = true;
	for (int i = 0; i <= M; i++) probv[i] = theta[i];
	for (int i = 0; i < nThreads; i++) {
		rc = pthread_create(&threads[i], &attr, E_STEP<ReadType, HitType, ModelType>, (void*)(&fparams[i]));
		pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) when calculating expected weights!");
	}
	for (int i = 0; i < nThreads; i++) {
		rc = pthread_join(threads[i], NULL);
		pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) when calculating expected weights!");
	}
	model.setNeedCalcConPrb(false);
	for (int i = 1; i < nThreads; i++) {
		for (int j = 0; j <= M; j++) {
			countvs[0][j] += countvs[i][j];
		}
	}
	countvs[0][0] += N0;

	/* destroy attribute */
	pthread_attr_destroy(&attr);


	sprintf(thetaF, "%s.theta", statName);
	fo = fopen(thetaF, "w");
	fprintf(fo, "%d\n", M + 1);

	// output theta'
	for (int i = 0; i < M; i++) fprintf(fo, "%.15g ", theta[i]);
	fprintf(fo, "%.15g\n", theta[M]);
	
	//calculate expected effective lengths for each isoform
	calcExpectedEffectiveLengths<ModelType>(model);
	polishTheta(theta, eel, model.getMW());

	// output theta
	for (int i = 0; i < M; i++) fprintf(fo, "%.15g ", theta[i]);
	fprintf(fo, "%.15g\n", theta[M]);

	fclose(fo);

	writeResults<ModelType>(model, countvs[0]);

	if (genBamF) {
		sprintf(outBamF, "%s.transcript.bam", outName);
		
		if (bamSampling) {
			READ_INT_TYPE local_N;
			HIT_INT_TYPE fr, to, len, id;
			vector<double> arr;
			uniform01 rg(engine_type(time(NULL)));

			if (verbose) cout<< "Begin to sample reads from their posteriors."<< endl;
			for (int i = 0; i < nThreads; i++) {
				local_N = hitvs[i]->getN();
				for (READ_INT_TYPE j = 0; j < local_N; j++) {
					fr = hitvs[i]->getSAt(j);
					to = hitvs[i]->getSAt(j + 1);
					len = to - fr + 1;
					arr.assign(len, 0);
					arr[0] = ncpvs[i][j];
					for (HIT_INT_TYPE k = fr; k < to; k++) arr[k - fr + 1] = arr[k - fr] + hitvs[i]->getHitAt(k).getConPrb();
					id = (arr[len - 1] < EPSILON ? -1 : sample(rg, arr, len)); // if all entries in arr are 0, let id be -1
					for (HIT_INT_TYPE k = fr; k < to; k++) hitvs[i]->getHitAt(k).setConPrb(k - fr + 1 == id ? 1.0 : 0.0);
				}
			}

			if (verbose) cout<< "Sampling is finished."<< endl;
		}

		BamWriter writer(inpSamType, inpSamF, pt_fn_list, outBamF, transcripts);
		HitWrapper<HitType> wrapper(nThreads, hitvs);
		writer.work(wrapper);
	}

	release<ReadType, HitType, ModelType>(readers, hitvs, ncpvs, mhps);
}

int main(int argc, char* argv[]) {
	ifstream fin;
	bool quiet = false;

	if (argc < 6) {
		printf("Usage : rsem-run-em refName read_type sampleName imdName statName [-p #Threads] [-b samInpType samInpF has_fn_list_? [fn_list]] [-q] [--gibbs-out] [--sampling]\n\n");
		printf("  refName: reference name\n");
		printf("  read_type: 0 single read without quality score; 1 single read with quality score; 2 paired-end read without quality score; 3 paired-end read with quality score.\n");
		printf("  sampleName: sample's name, including the path\n");
		printf("  sampleToken: sampleName excludes the path\n");
		printf("  -p: number of threads which user wants to use. (default: 1)\n");
		printf("  -b: produce bam format output file. (default: off)\n");
		printf("  -q: set it quiet\n");
		printf("  --gibbs-out: generate output file used by Gibbs sampler. (default: off)\n");
		printf("  --sampling: sample each read from its posterior distribution when bam file is generated. (default: off)\n");
		printf("// model parameters should be in imdName.mparams.\n");
		exit(-1);
	}

	time_t a = time(NULL);

	strcpy(refName, argv[1]);
	read_type = atoi(argv[2]);
	strcpy(outName, argv[3]);
	strcpy(imdName, argv[4]);
	strcpy(statName, argv[5]);

	nThreads = 1;

	genBamF = false;
	bamSampling = false;
	genGibbsOut = false;
	pt_fn_list = pt_chr_list = NULL;

	for (int i = 6; i < argc; i++) {
		if (!strcmp(argv[i], "-p")) { nThreads = atoi(argv[i + 1]); }
		if (!strcmp(argv[i], "-b")) {
			genBamF = true;
			inpSamType = argv[i + 1][0];
			strcpy(inpSamF, argv[i + 2]);
			if (atoi(argv[i + 3]) == 1) {
				strcpy(fn_list, argv[i + 4]);
				pt_fn_list = (char*)(&fn_list);
			}
		}
		if (!strcmp(argv[i], "-q")) { quiet = true; }
		if (!strcmp(argv[i], "--gibbs-out")) { genGibbsOut = true; }
		if (!strcmp(argv[i], "--sampling")) { bamSampling = true; }
	}

	general_assert(nThreads > 0, "Number of threads should be bigger than 0!");

	verbose = !quiet;

	//basic info loading
	sprintf(refF, "%s.seq", refName);
	refs.loadRefs(refF);
	M = refs.getM();
	sprintf(groupF, "%s.grp", refName);
	gi.load(groupF);
	m = gi.getm();

	sprintf(tiF, "%s.ti", refName);
	transcripts.readFrom(tiF);

	sprintf(cntF, "%s.cnt", statName);
	fin.open(cntF);

	general_assert(fin.is_open(), "Cannot open " + cstrtos(cntF) + "! It may not exist.");

	fin>>N0>>N1>>N2>>N_tot;
	fin.close();

	general_assert(N1 > 0, "There are no alignable reads!");

	if ((READ_INT_TYPE)nThreads > N1) nThreads = N1;

	//set model parameters
	mparams.M = M;
	mparams.N[0] = N0; mparams.N[1] = N1; mparams.N[2] = N2;
	mparams.refs = &refs;

	sprintf(mparamsF, "%s.mparams", imdName);
	fin.open(mparamsF);

	general_assert(fin.is_open(), "Cannot open " + cstrtos(mparamsF) + "It may not exist.");

	fin>> mparams.minL>> mparams.maxL>> mparams.probF;
	int val; // 0 or 1 , for estRSPD
	fin>>val;
	mparams.estRSPD = (val != 0);
	fin>> mparams.B>> mparams.mate_minL>> mparams.mate_maxL>> mparams.mean>> mparams.sd;
	fin>> mparams.seedLen;
	fin.close();

	//run EM
	switch(read_type) {
	case 0 : EM<SingleRead, SingleHit, SingleModel>(); break;
	case 1 : EM<SingleReadQ, SingleHit, SingleQModel>(); break;
	case 2 : EM<PairedEndRead, PairedEndHit, PairedEndModel>(); break;
	case 3 : EM<PairedEndReadQ, PairedEndHit, PairedEndQModel>(); break;
	default : fprintf(stderr, "Unknown Read Type!\n"); exit(-1);
	}

	time_t b = time(NULL);

	printTimeUsed(a, b, "EM.cpp");

	return 0;
}
