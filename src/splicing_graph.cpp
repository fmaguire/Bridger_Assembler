#include "splicing_graph.h"
#include <iostream>
#include <algorithm>

bool is_similar(const std::string& str1, const std::string& str2, char mode) {

  int mismatch = 0;
  // note , str1 and str2 should have corresponding bases from begein or end
  int length = str1.length() > str2.length() ? str2.length(): str1.length();

  if (length == 0) { 
    if (str1.empty() && str2.empty())
      return true;
    else
      return false;
  }

  if (mode == 'F') {
    for (int i = 0; i < length; ++i) {
      if (str1[i] != str2[i])
        mismatch++; 
    }
  } else {
    for (int i = 0; i < length; ++i) {
      if (str1[str1.length()-i-1] != str2[str2.length()-i-1])
        mismatch++;
    }
  }

  if ((float)mismatch/length < 0.35) {
    return true;
  } else {
    return false;
  }
}


bool is_aligned(const std::string& str1, const std::string& str2) {

  int mismatch = 0;

  // note , str1 and str2 should have corresponding bases from begein
  int length = str1.length() > str2.length() ? str2.length(): str1.length();
  
  for (int i = 0; i < length; ++i) {
    if (str1[i] != str2[i]) 
      mismatch++; 
  }

  return (mismatch <= 2);
}


bool compatible(const std::string& str1, const std::string& str2) {  
  
  for (unsigned int i = 0; i <= str2.length()-g_kmer_length; ++i) {
    const std::string& kmer = str2.substr(i,g_kmer_length);
    std::string::size_type start = str1.find(kmer);
    if (start != std::string::npos) {
      if (start > i)
	return is_aligned(str1.substr(start-i),str2);
      else 
	return is_aligned(str1, str2.substr(i-start));
    }
  }

  return false;
}

/*

std::string get_consensus(const std::string& kmer, Kmer_hash_map& kmer_hash,
                           std::vector<std::string>& reads, int& breakpoint){

    breakpoint = -1;
    int flag;   // error type : mismatch-0 , deletion-1, insertion-2

    std::string temp = reads[0];
    std::vector<int> cov;

cout << "temp=" << temp << endl;

    for (int i= 0; i < temp.size(); ++i) {
        cov.push_back(1);
    }

    for (unsigned int j = 1; j< reads.size(); ++j) {

        std::string str = reads[j];
cout << "str =" << str << endl;

  int length = str.size() < temp.size() ? str.size() : temp.size();
        for (int k = 0; k < length; ++k) {

      if (str[k] == temp[k]) {
      cov[k]++;
      }
      else {

    int dist = length - k -1;
cout << "dist=" << dist << endl;

    if (dist < 1 || (str.size()-k) < 5) {
      break;
    }
    else if(dist >10){
        dist = 10;  // 10 base is enough
    }

    if (temp.substr(k+1,dist) == str.substr(k+1,dist)){
        flag = 0;    // a mismatch
        correction(temp, str, kmer_hash, kmer, k, flag);
        continue;
    }

    std::string::size_type idx1 = temp.substr(k,dist).find(str.substr(k+1,dist));
    std::string::size_type idx2 = str.substr(k,dist).find(temp.substr(k+1,dist));
    if( idx1 == std::string::npos){
        flag = 1;   // deletion
        correction(temp, str, kmer_hash, kmer, k, flag);
        break;
    }
    else if( idx2 == std::string::npos) {
        flag = 2;      // insertion
        correction(temp, str, kmer_hash, kmer, k, flag);
        break;
    }
    else{
        breakpoint = k;
    }
            }
  }

  while(length < str.length()) {
      temp += str.substr(length-1,1);
      cov.push_back(1);
      ++length;
  }
cout << "temp=" << temp << endl;
    }

    //
    int k = 0;
    for(; k < temp.size(); ++k) {
  cerr <<cov[k]<< " ";
  if(cov[k] < g_min_kmer_coverage){
      break;
  }
    }

    cerr << endl;
    if (breakpoint > 0 && breakpoint < k){
  return temp.substr(0,breakpoint+1);
    }
    return temp.substr(0,k+1);
}


bool correction(std::string& s1, std::string& s2, Kmer_hash_map& kmer_hash,
      const std::string kmer, int pos, int flag) {

    std::string str1 = kmer + s1;
    std::string str2 = kmer + s2;
    std::string kmer1 = str1.substr(pos+1, g_kmer_length);
    std::string kmer2 = str2.substr(pos+1, g_kmer_length);
    kmer_int_type_t intval1 = kmer_to_intval(kmer1);
    kmer_int_type_t intval2 = kmer_to_intval(kmer2);

    if(kmer_hash[intval1].size() < kmer_hash[intval2].size()) {
  if (flag == 0) {
      s1[pos] = s2[pos];
  }
  else if(flag == 1) {
      std::string temp = s1.substr(0,pos) + s2[pos] + s1.substr(pos);
      s1 = temp;
  }
  else{
      s1.erase(pos);
  }
  return true;
    }
    return false;
}

*/

