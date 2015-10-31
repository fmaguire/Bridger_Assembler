#include "loadreads.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>

const int MAX_STR = 1024;

void load_data(std::vector<std::string>& data, std::string file) {

  data.reserve(1000000);
  std::fstream in;
  in.open(file.c_str(), std::fstream::in);
  if (!in.is_open()) {
    std::cerr << "[error] File " << file << " can't be opened." << std::endl;
    exit(1);
  }

  std::cerr << "Begin loading reads ..." << std::endl;

  char c_line[MAX_STR];
  while(!in.eof()) {
    in.getline(c_line, MAX_STR);
    if (c_line[0] == '>') {
      in.getline(c_line, MAX_STR);
      std::string sequence(c_line);
      data.push_back(sequence);
    }
  }

  std::cerr << data.size() <<" reads have been loaded !" << std::endl;
  in.close();
}


