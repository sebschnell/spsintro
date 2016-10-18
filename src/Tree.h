#ifndef TREE_H // header guard to copy the content of the .h file only once into the program
#define TREE_H
#include <string>

struct Tree {
  std::string spec_name;
  unsigned int no_stems = 0;
  double height = 0.0;
  double dbh = 0.0;
};
#endif