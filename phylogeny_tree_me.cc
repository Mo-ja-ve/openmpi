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

void broadcast(vector <char> &char_string);

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

  int num_genomes = genomes.size();
  MPI_Bcast(&num_genomes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  genomes.resize(num_genomes);

     for(int i = 0; i < genomes.size(); i++){

       vector <char> char_string(genomes[i].begin(), genomes[i].end());
       broadcast(char_string);
       string temp(char_string.begin(),char_string.end());
       genomes[i] = temp;
     }

  cout <<endl<< "Number of covid genomes = " << genomes.size()<< endl;

  // Make initial labels on all strings
  vector<pair<string,string> > genome_tree;
  for (int i=0;i<genomes.size();i++) {
    genome_tree.push_back(make_pair(to_string(i),genomes[i]));
  }


  while (genome_tree.size() >1) {

       vector <pair<int,int>> proc_pair;

       for(int i = 0; i <genome_tree.size(); i++){
            for(int j = i+1; j < genome_tree.size(); j++){
                 proc_pair.push_back(make_pair(i,j));
            }
       }
  int max_i = 0;
  int max_j = 0;
  string best;
  vector <char> char_best;

  int proc_longestLCS;
  int proc_smallest_i;
  bool start = true;

  int loc_length, global_longest;
  if(myid !=0){
    for(int k = myid-1; k < proc_pair.size(); k = k + num_procs-1){
         //cout<<endl<<"K: "<<k;
      int i = proc_pair[k].first;
      int j = proc_pair[k].second;

      string z;
      z=compute_LCS(genome_tree[i].second, genome_tree[j].second);
      loc_length=z.length();
      max_i = i;
      max_j = j;
      MPI_Send(&loc_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&max_i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&max_j, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);


     //MPI_Allreduce(&j, &max_j, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

           //max_j = j;
           // cout<<endl<<"MAX I: "<<max_i;
           //cout<<endl<<"I: "<<i;
           //cout<<endl<<"J: "<<j;
           // cout<<endl<<"hello! ";
           // cout<<endl<<"MAX J: "<<max_j;
           // copy(z.begin(), z.end(), back_inserter(char_best));
           // int size2[2];
           // size2[0] = char_best.size();

           // MPI_Bcast(size2, 1, MPI_INT, myid, MPI_COMM_WORLD);
           // char_best.resize(size2[0]);
           // MPI_Bcast(&char_best[0], size2[0], MPI_CHAR, myid, MPI_COMM_WORLD);
           //for (const char &c: char_best)
               //std::cout << c;
          // string best_temp(char_best.begin(), char_best.end());
          // best = best_temp;
     }

     }
     }


     if(myid == 0){

          vector <vector <int>> largest_lcs;
          largest_lcs.resize(num_procs-1);
          for(int j = 0; j < largest_lcs.size(); j++){
               largest_lcs[j].resize(3);
          }

          for(int i = 1; i < num_procs; i++){
               MPI_Recv(&loc_length, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&max_i, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               MPI_Recv(&max_j, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

               largest_lcs[i-1][0] = loc_length;
               largest_lcs[i-1][1] = max_i;
               largest_lcs[i-1][2] = max_j;

          }
          int temp_lcs;
          int temp_max_i;
          int temp_max_j;

          for(int i = 0; i <largest_lcs.size(); i++){
               for(int j = i+1; j < largest_lcs.size(); j++){
                    if(largest_lcs[i][0] < largest_lcs[j][0]){
                         temp_lcs   = largest_lcs[i][0];
                         temp_max_i = largest_lcs[i][1];
                         temp_max_j = largest_lcs[i][2];

                         largest_lcs[i][0] = largest_lcs[j][0];
                         largest_lcs[i][1] = largest_lcs[j][1];
                         largest_lcs[i][2] = largest_lcs[j][2];

                         largest_lcs[j][0] = temp_lcs;
                         largest_lcs[j][1] = temp_max_i;
                         largest_lcs[j][2] = temp_max_j;

                    }
               }
          }

          max_i = largest_lcs[0][1];
          max_j = largest_lcs[0][2];
          // cout<<endl<<"max i: "<<max_i;
          // cout<<endl<<"max j: "<<max_j;

          for(int i = 1; i < num_procs; i++){
               MPI_Send(&max_i, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
               MPI_Send(&max_j, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
          }
     }

     if(myid !=0){
     MPI_Recv(&max_i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     MPI_Recv(&max_j, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

     // cout<<endl;
     // cout<<endl<<"my proc id"<<myid<<"max i: "<<max_i;
     // cout<<endl<<"my proc id"<<myid<<"max j: "<<max_j;

     best = compute_LCS(genome_tree[max_i].second, genome_tree[max_j].second);

     string new_tree_label = "("+genome_tree[max_i].first + "," + genome_tree[max_j].first +")";
     genome_tree.erase(genome_tree.begin()+max_i);
     genome_tree.erase(genome_tree.begin()+max_j-1); // max_i got deleted!
     genome_tree.push_back(make_pair(new_tree_label,best));
     cout<<endl<<"size :"<<genome_tree.size();
     }


     if(myid == 1){
          cout << "Phylogeny = " << endl;
          cout << genome_tree[0].first << endl;
          cout << "Root has length " << genome_tree[0].second.length() << endl;
          // cout << "Best pair = " << max_i << " " << max_j << endl;
          // cout << "Length = " << best.length() << endl;
     }

  // Debug
  // for (int i=0;i<genomes.size();i++) {
  //   assert(Is_subsequence(0,0,genome_tree[0].second,genomes[i]));
  // }

 MPI_Finalize();
}

void broadcast(vector <char> &char_string){

  //  set all this up based off slide notes
  int buffer[2];
  buffer[0] = char_string.size();

  MPI_Bcast(buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);

  char_string.resize(buffer[0]);

   MPI_Bcast(&char_string[0], buffer[0], MPI_CHAR, 0, MPI_COMM_WORLD );
}
