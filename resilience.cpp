#include "common.h"
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <cmath>

using namespace std;

// Helper that analyzes fragmentation for a single node removal strategy
// Calculates distribution of strongly connected component sizes
    int res(int organism, string strategy,bool simple_output) {
    //static const int NUM_POINTS = 100;
	int NUM_POINTS = 100;
    string input_path = get_output_path(FRAGMENT_DIR, organism, simple_output);
    input_path += strategy;
    ifstream fragmentation(input_path.c_str(), ios_base::app);
	
	string output_path = get_output_path(FRAGMENT_DIR, organism, simple_output);
    output_path += strategy;
	output_path += "hmsh";
	ofstream hmshfile(output_path.c_str());
	
	double f=0;
	int nodesum=0;
	int fragNode=0;
	int fragNum=0;
	string line;
	int index=0;
	string y;
	double pi=0;
	double pi2=0;
	double pilog=0;
  	double hmsh=0;
	double m=0;
	double res=0;
	string x ;
	double pilog2=0;
	
	fragmentation.clear();
	fragmentation.seekg(0);
	
	hmshfile << fixed << setprecision(PRECISION);
	

	getline(fragmentation, line);
	sscanf(line.c_str(), "%d", &NUM_POINTS); 
	cout<<"NUM_POINTS:"<<NUM_POINTS<<endl; 
	
	while (getline(fragmentation, line)) {
		pi=0;
		pi2=0;
		pilog=0;
		pilog2=0;
		
		index = line.find(" ",0);
		if(index < line.length())
			;//cout<<"Found at index : "<< index <<endl;
		else
			;//cout<<"Not found"<<endl;
		x = line.substr(0,index);
		//cout<<x<<endl; 
	
		size_t size = line.size();
		
		y= line.substr(index+1,size);
		
		sscanf(y.c_str(), "%d", &nodesum);
		
		if(nodesum>=1)
		{
			if(nodesum>1) 
			{
				index = y.find(",",0);
				size = y.size();

				while((size-index>2) &&(index >0))
				{ 
					y= y.substr(index+1,size);
					sscanf(y.c_str(), "%d %d", &fragNode, &fragNum);
					pi =(double)fragNode/NUM_POINTS;
					pilog +=pi*log(pi)*fragNum;

					
					index = y.find(",",0);
					size = y.size();		

				}
			}

		if(nodesum == 1) 
			nodesum=0;
			
		pi2=(double)(NUM_POINTS-nodesum)/NUM_POINTS;
		
		pilog2 = pi2*log((double)1/NUM_POINTS);
		
		pilog +=pilog2;
		
		m=-1*pilog/log(NUM_POINTS);   //node_removal 100

		
		hmsh += m;  

		hmshfile<<x<<" "<<m<<endl;
		
		}
		
	}
	res = 1-hmsh/100; 
	cout<<"res:"<<res<<endl; 
	return 0;
}

// Creates a graph out of the organism's data, where vertices are proteins and edges are interactions
PGraph create_graph(int organism, string input, vector<string> &proteins) {
    string data_path = (input == "") ? get_path(EXTRACTED_DATA_DIR, organism) : input;
    ifstream data(data_path.c_str());
    string line;
    char first_protein[MAX_LENGTH];
    char second_protein[MAX_LENGTH];
	
    PGraph G = PGraph::TObj::New();
    // We need a map in order to retrieve the correct node for a protein
    unordered_map<string, TInt> protein_to_node;


    while (getline(data, line)) {
        sscanf(line.c_str(), "%[^ ] %s", first_protein, second_protein);
        if (protein_to_node.count(first_protein) == 0) {
            protein_to_node[first_protein] = G->AddNode();
            proteins.push_back(first_protein);
        }
        if (protein_to_node.count(second_protein) == 0) {
            protein_to_node[second_protein] = G->AddNode();
            proteins.push_back(second_protein);
        }
        G->AddEdge(protein_to_node[first_protein], protein_to_node[second_protein]);


    }

    IAssert(G->IsOk());
    return G;
}

void EdgeRemove(int organism, string input,int removeRate) {
	vector<string> proteins;
	TIntFltH NIdRndH;
	FILE * pfile1=fopen("./EdgeRemove_before","wb");
	FILE * pfile2=fopen("./EdgeRemove_later","wb");
	
	PGraph G = create_graph(organism, input, proteins); 
	
	TRnd random(0); 
	TRnd Rnd=1;
	
	G = TSnap::GetMxScc(G); 	
	TSnap::SaveEdgeList(G, "graphold.txt", "Save as tab-separated list of edges");

	G->Dump(pfile1);
	fclose(pfile1);
	
	cout<<"nodes:"<<G->GetNodes()<<" edges:"<<G->GetEdges()<<endl;
	
	PGraph H = new PGraph::TObj(*G);


	H = TSnap::GenRewire(G, removeRate);
	cout<<"nodes:"<<H->GetNodes()<<", edges:"<<H->GetEdges()<<endl;
	
	TSnap::SaveEdgeList(H, "graphnew1.txt", "Save as tab-separated list of edges");
	H->Dump(pfile2);
	fclose(pfile2);
	

}

int main(int argc, char* argv[]) {
    int organism = atoi(argv[1]);

	res(organism, argv[2], false); // remove random

    return 0;
}