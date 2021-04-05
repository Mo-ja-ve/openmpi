// A more memory efficient version of phylogeny_tree.cc
// Early version consumed 85% of RAM on a 16GB machine.
//
//
// A Basic C++ implementation of
// a program that computes a Phylogenetic tree from COVID-19 Virus sequences
// Each virus sequence has about 30K base pairs.
// Computing the LCS takes about 7 seconds.
//
// Basic algorithm to build the phylogeny:
// Compute all LCS (longest common subsequence) for all pairs of
// strings.
// Pick the best pair.
// Merge the best into a new string that is the LCS of those two.
// Continue until only one sequence remains.
//
// Written by David W. Juedes
//
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <mpi.h>
#include <fstream>

using namespace std;

void broadcast(vector <const char*> char_string);

// This code checks whether
// x[i..n-1] a sequence of y[j..m-1], where n is the length of
// x, m is the length of y.
// This is solely used to check the correctness of the solution.

bool Is_subsequence(int i, int j, string &x, string &y) {
  if (i>=x.length()) return true; // Consumed all of x's characters.
  if (j>=y.length()) return false; // Can't consume y's first.

  if (x[i]==y[j]) {
    return Is_subsequence(i+1,j+1,x,y);
  } else {
    return Is_subsequence(i,j+1,x,y);
  }
}


int max3(int x,int y, int z) {
  return max(max(x,y),z);
}

int max3_loc(int x,int y, int z) {
  if (max3(x,y,z) == x) return 1;
  if (max3(x,y,z) == y) return 2;
  if (max3(x,y,z) == z) return 3;
  return -1;
}
int max_loc(int x,int y) {
  if (max(x,y) == x) return 1;
  else {return 2;}
}

//
// Recursively construct the longest common subsequence from the dynamic
// programming table "from"
//
string rec_string(string &x1, string &y1, vector<vector<int> > &LCS,
		  //vector<vector<pair<int,pair<int,int> > > > & from,
		  int m,
		  int n) {

  //  pair<int,pair<int,int> > t;
  if (n==0) return "";
  if (m==0) return "";
  int a = LCS[m-1][n-1];
  int b = LCS[m-1][n];
  int c = LCS[m][n-1];
  // Remember, we are off by 1.
  if (x1[m-1] == y1[n-1]) {
    if (((a+1) >=b) && ((a+1) >=c)) {
      return rec_string(x1,y1,LCS,m-1,n-1)+x1[m-1];
    }
    else {
      if (b>=c) {
	return rec_string(x1,y1,LCS,m-1,n);
      } else {
	return rec_string(x1,y1,LCS,m,n-1);
      }
    }
  } else {
          if (b>=c) {
	    return rec_string(x1,y1,LCS,m-1,n);
	  } else {
	    return rec_string(x1,y1,LCS,m,n-1);
	  }
  }

}


string compute_LCS(string &x1, string &y1) {
  vector<vector<int> > LCS;
  // vector<vector<pair<int, pair<int,int> > > > from;

  //from.resize(x1.length()+1);
  LCS.resize(x1.length()+1);
  for (int i=0;i<=x1.length();i++) {
    LCS[i].resize(y1.length()+1);
    //from[i].resize(y1.length()+1);
  }
  LCS[0][0] = 0;
  //from[0][0] = make_pair(0,make_pair(-1,-1));
  for (int i=0;i<=x1.length();i++) {
    LCS[i][0]=0;
    //from[i][0] = make_pair(0,make_pair(-1,-1));
  }
  for (int i=0;i<=y1.length();i++) {
    LCS[0][i]=0;
    //from[0][i] = make_pair(0,make_pair(-1,-1));
  }

  for (int i=1;i<=x1.length();i++) {
    for (int j=1;j<=y1.length();j++) {
      if (x1[i-1]==y1[j-1]) {

	LCS[i][j] = max3(LCS[i-1][j-1]+1,LCS[i-1][j],LCS[i][j-1]);
	/*	switch (max3_loc(LCS[i-1][j-1]+1,LCS[i-1][j],LCS[i][j-1])) {
	case 1: from[i][j] = make_pair(1,make_pair(i-1,j-1));
	  break;
	case 2: from[i][j] = make_pair(0,make_pair(i-1,j));
	  break;
	case 3: from[i][j] = make_pair(0,make_pair(i,j-1));
	  break;
	  } */
      } else {
	LCS[i][j] = max(LCS[i-1][j],LCS[i][j-1]);
	/* switch (max_loc(LCS[i-1][j],LCS[i][j-1])) {
	case 1: from[i][j] = make_pair(0,make_pair(i-1,j));
	  break;
	case 2: from[i][j] = make_pair(0,make_pair(i,j-1));
	  break;
	  } */
	}
    }
  }
  //cout << LCS[x1.length()][y1.length()] << endl;
  string z = rec_string(x1,y1,LCS,x1.length(),y1.length());

  assert(z.length() == LCS[x1.length()][y1.length()]);
  assert(Is_subsequence(0,0,z,x1));
  assert(Is_subsequence(0,0,z,y1));

  return z;
}




/// Compute the longest common subsequence.
int main(int argc, char *argv[]) {

  fstream fin;

  //  mpi init stuff
  int num_procs, myid;
  double timer_total;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  //cout<<endl<<"num proc: "<<num_procs;
  //cout<<endl<<"my id: "<<myid;

  vector<string> genomes;
  if(myid == 0){

    while (!cin.eof()) {
      string line;
      getline(cin,line);
      if (!cin.fail()) {
        genomes.push_back(line);
      }
    }
  }

  vector <const char*> char_string(genomes.size(), nullptr);

  for(int i = 0; genomes.size(); i++){
       char_string[i] = genomes[i].c_str();
 }

  broadcast(char_string);

  cout <<endl<< "Number of covid genomes = " << genomes.size()<< endl;

  // Make initial labels on all strings
  vector<pair<string,string> > genome_tree;
  for (int i=0;i<genomes.size();i++) {
    genome_tree.push_back(make_pair(to_string(i),genomes[i]));
  }
  cout<<endl<<"MY ID: "<<myid;
  //
  //
  //
  // while (genome_tree.size() >1) {
  //      vector <pair<int,int>> proc_pair;
  //
  //      for(int i = 0; i <genomes.size()-1; i++){
  //           for(int j = i+1; j < genomes.size(); j++){
  //                proc_pair.push_back(make_pair(i,j));
  //           }
  //      }
  //
  //      pair<int,int> recv_ij;
  //      for(int i = myid; i < proc_pair.size(); i+=num_procs){
  //
  //           MPI_Allreduce(&proc_pair[i], &recv_ij, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  //      }
  //
  //     cout<<endl<<"RECV_IJ: "<<recv_ij.first<<" "<<recv_ij.second;
  //
  //   // Debugging
  //   //cout << "Tree size = " << genome_tree.size() << endl;
  //
  // int max_i;
  // int max_j;
  // string best;
  //
  // bool start = true;
  // bool start_recv = start;
  // int loc_length, recv_length;
  //
  //   //cout << "i = " << i << endl;
  //     string z;
  //     z=compute_LCS(genome_tree[recv_ij.first].second, genome_tree[recv_ij.second].second);
  //
  //     if (start_recv) {
	// start = false;
  //    MPI_Allreduce(&start, &start_recv, 1 , MPI_C_BOOL, MPI_MIN, MPI_COMM_WORLD);
	// max_i = recv_ij.first;
	// max_j = recv_ij.second;
	// best = z;
  //     } else {
  //
  //    loc_length = z.length();
  //    MPI_Allreduce(&loc_length, &recv_length, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  //
  //    if (recv_length > best.length()) {
	//   max_i = recv_ij.first;
	//   max_j= recv_ij.second;
	//   best = z;
	// }
  // }
  // string new_tree_label = "("+genome_tree[max_i].first + "," + genome_tree[max_j].first +")";
  // genome_tree.erase(genome_tree.begin()+max_i);
  // genome_tree.erase(genome_tree.begin()+max_j-1); // max_i got deleted!
  // genome_tree.push_back(make_pair(new_tree_label,best));
  // }



// cout << "Phylogeny = " << endl;
// cout << genome_tree[0].first << endl;
// cout << "Root has length " << genome_tree[0].second.length() << endl;
// cout << "Best pair = " << max_i << " " << max_j << endl;
// cout << "Length = " << best.length() << endl;


  // Debug
  // for (int i=0;i<genomes.size();i++) {
  //   assert(Is_subsequence(0,0,genome_tree[0].second,genomes[i]));
  // }

 MPI_Finalize();
}

void broadcast(vector <const char*> char_string){

  //  set all this up based off slide notes
  int buffer[2];
  buffer[0] = char_string.size();

  MPI_Bcast(buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);

  char_string.resize(buffer[0]);

   MPI_Bcast(&char_string[0], buffer[0], MPI_CHAR, 0, MPI_COMM_WORLD );
}
