#include "common.h"
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <cmath>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
double res(int organism, string strategy,bool simple_output) ;
int beforeAnalyzeRemove(PGraph G,int removepercent);
void NodeSortByRandom(PGraph G);

TIntFltH NIdRndH11;  
TIntFltH::TIter TISort;

// Creates the necessary directories
void make_analyze_dir() {
	cout << make_analyze_dir  << endl;
    make_dir_rooted(STAT_DIR);
    make_dir_rooted(KCORE_DIR);
    make_dir_rooted(DEGREE_DIR);
    make_dir_rooted(CURVE_DIR);
    make_dir_rooted(CLUSTER_DIR);
    make_dir_rooted(LORENZ_DIR);
    make_dir_rooted(DEGREECTR_DIR);
    make_dir_rooted(BETWEEN_DIR);
    make_dir_rooted(CLOSE_DIR);
    //make_dir_rooted(ADJACENCY_DIR);
    make_dir_rooted(FRAGMENT_DIR);
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
		//cout<<first_protein<<", "<<second_protein<<endl;
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

// Calculates average degree
double average_degree(PGraph G) {
    int total_degree = 0;
    int nodes = G->GetNodes();
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        total_degree += NI.GetDeg();
    }
    return (double)total_degree / nodes;
}

// Calculates maximum degree
int maximum_degree(PGraph G) {
    int max_degree = 0;
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        int curr_degree = NI.GetDeg();
        if (curr_degree > max_degree) max_degree = curr_degree;
    }
    return max_degree;
}

// Calculates density, (number of edges present) / (number of edges possible)
double density(PGraph G) {
    int edges_present = G->GetEdges();
    int nodes = G->GetNodes();
    int edges_possible = nodes * (nodes - 1) / 2;
    return (double)edges_present / edges_possible;
}

// Calculates number of connected components
int num_components(PGraph G) {
    TCnComV CnComV;
    TSnap::GetSccs(G, CnComV);
    return CnComV.Len();
}

// Calculates global clustering coefficient, (number of closed triads) / (number of triads)
double global_clustering(PGraph G) {
    int64 closed = 0, open = 0;
    TSnap::GetTriads(G, closed, open);
    return (double)closed / (closed + open);
}

// Returns n choose k
int64 choose(int n, int k) {
    int64 result = 1;
    for (int i = 0; i < k; i++) {
        // Separate into two lines to prevent truncation during integer division
        result *= n - i;
        result /= (i + 1);
    }
    return result;
}

// Calculates density of k-stars, (number of k-stars present) / (number of k-stars possible)
double star_density(PGraph G, int k) {
    int64 count = 0;
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        count += choose(NI.GetDeg(), k);
    }
    int nodes = G->GetNodes();
    // Divide in two steps to prevent overflow
    double density = (double)count / nodes;
    return density / choose(nodes, k);
}

// Calculates the Gini coefficient, a measure of degree inequality
double gini_coefficient(PGraph G) {
    TIntPrV DegToCntV;
    TSnap::GetDegCnt(G, DegToCntV);

    int nodes = G->GetNodes();
    int64 moment_degree_sum = 0; // sum of i * d_i over a sorted list of degrees d_1, d_2, ..., d_n
    int64 degree_sum = 2 * G->GetEdges();

    int nodes_so_far = 0;

    for (TIntPrV::TIter TI = DegToCntV.BegI(); TI < DegToCntV.EndI(); TI++) {
        TInt degree = TI->GetVal1();
        TInt frequency = TI->GetVal2();

        // We want to sum the series (nodes_so_far + 1) + (nodes_so_far + 2) + ... + (nodes_so_far + frequency)
        int64 series_sum = frequency * nodes_so_far + frequency * (frequency + 1) / 2;

        moment_degree_sum += series_sum * degree;
        nodes_so_far += frequency;
    }

    double result = 2 * moment_degree_sum / ((double)nodes * degree_sum);
    result -= ((double)nodes + 1) / nodes;
    return result;
}

// Calculates the normalized edge distribution entropy, a measure of degree equality
double edge_entropy(PGraph G) {
    TIntPrV DegToCntV;
    TSnap::GetDegCnt(G, DegToCntV);

    int nodes = G->GetNodes();
    int edges = G->GetEdges();
    double result = 0;

    for (TIntPrV::TIter TI = DegToCntV.BegI(); TI < DegToCntV.EndI(); TI++) {
        TInt degree = TI->GetVal1();
        TInt frequency = TI->GetVal2();
        double degree_fraction = (double)degree / (2 * edges);
        result -= frequency * degree_fraction * log(degree_fraction);
    }

    result /= log(nodes);
    return result;
}

// Calculates the Pearson correlation of the degrees at the ends of edges, a measure of assortative mixing
double assortative_mixing(PGraph G) {
    double inverse_of_edges = (double)1 / G->GetEdges();
    long long sum_of_products = 0, sum_of_sums = 0, sum_of_squares = 0;
    for (PGraph::TObj::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
        int src = G->GetNI(EI.GetSrcNId()).GetDeg();
        int dst = G->GetNI(EI.GetDstNId()).GetDeg();
        sum_of_products += src * dst;
        sum_of_sums += src + dst;
        sum_of_squares += src * src + dst * dst;
    }

    double sum_of_products_term = inverse_of_edges * sum_of_products;
    double sum_of_sums_term = inverse_of_edges * sum_of_sums / 2;
    sum_of_sums_term *= sum_of_sums_term;
    double sum_of_squares_term = inverse_of_edges * sum_of_squares / 2;
    return (sum_of_products_term - sum_of_sums_term) / (sum_of_squares_term - sum_of_sums_term);
}

// Analyzes basic summary statistics of the graph
void analyze_statistics(PGraph G, int organism, bool simple_output) {
    string output_path = get_output_path(STAT_DIR, organism, simple_output);
    //ofstream statistics(output_path.c_str(), ios_base::app);
    ofstream statistics(output_path.c_str());

    statistics << fixed << setprecision(PRECISION)
        << G->GetNodes() << endl
        << G->GetEdges() << endl
        << average_degree(G) << endl
        << maximum_degree(G) << endl
        << density(G) << endl
        << num_components(G) << endl
        << TSnap::GetMxSccSz(G) << endl
        << TSnap::GetBfsFullDiam(G, 100) << endl
        << TSnap::GetBfsEffDiam(G, 100) << endl
        << global_clustering(G) << endl
        << TSnap::GetClustCf(G) << endl
        << star_density(G, 2) << endl
        << star_density(G, 3) << endl
        << gini_coefficient(G) << endl
        << edge_entropy(G) << endl
        << assortative_mixing(G) << endl;
}

// Analyzes the k-core distribution of the graph
void analyze_kcores(PGraph G, int organism, bool simple_output) {
    string output_path = get_output_path(KCORE_DIR, organism, simple_output);
    ofstream kcores(output_path.c_str());
    TKCore<PGraph> KCore(G);

    while (KCore.GetNextCore()) {
        kcores << KCore.GetCurK() << " "
            << KCore.GetCoreNodes() << " "
            << KCore.GetCoreEdges() << endl;
    }
}

//Analyzes node degree
void analyze_nodedegree(PGraph G, int organism, bool simple_output,const vector<string> &proteins) {
	string output_path = get_output_path("nodedegree", organism, simple_output);
	ofstream nodedegree(output_path.c_str()); 
	
	for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        int NId = NI.GetId();
        nodedegree << proteins.at(NId) << " "<< NId << " "<< NI.GetDeg() << endl;

    }
}
// Analyzes the degree distribution of the graph
void analyze_degrees(PGraph G, int organism, bool simple_output) {
    string output_path = get_output_path(DEGREE_DIR, organism, simple_output);
    ofstream degrees(output_path.c_str());
    TIntPrV DegToCntV;
	
    TSnap::GetDegCnt(G, DegToCntV); 

    for (TIntPrV::TIter TI = DegToCntV.BegI(); TI < DegToCntV.EndI(); TI++) {
        degrees << TI->GetVal1() << " "
            << TI->GetVal2() << endl;
    }
}

// Analyzes the node curvature distribution of the graph
void analyze_curvature(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins) {
    string output_path = get_output_path(CURVE_DIR, organism, simple_output);
    ofstream curvature(output_path.c_str());

    curvature << fixed << setprecision(PRECISION);

    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        int NId = NI.GetId();
        double node_curvature = 1 - (double)NI.GetDeg() / 2 + (double)TSnap::GetNodeTriads(G, NId) / 3;
        curvature << proteins.at(NId) << " "<< NId << " "<< node_curvature << endl;

    }
}

// Analyzes the clustering coefficient distribution of the graph
void analyze_clustering(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins) {
    string output_path = get_output_path(CLUSTER_DIR, organism, simple_output);
    ofstream clustering(output_path.c_str());
    TIntFltH NIdCCfH;
    TSnap::GetNodeClustCf(G, NIdCCfH);

    clustering << fixed << setprecision(PRECISION);

    for (TIntFltH::TIter TI = NIdCCfH.BegI(); TI < NIdCCfH.EndI(); TI++) {
        clustering << proteins.at(TI.GetKey()) << " "<< TI.GetKey() << " "<< TI.GetDat() << endl;
    }
}
//Analyzes Betweenness Center
void analyze_GetBetweennessCentr(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins) {
    string output_path = get_output_path("betweencenter", organism, simple_output);
    ofstream betweencenter(output_path.c_str());
    TIntFltH NIdBCfH;
    TSnap::GetBetweennessCentr(G, NIdBCfH);

    betweencenter << fixed << setprecision(PRECISION);

    for (TIntFltH::TIter TI = NIdBCfH.BegI(); TI < NIdBCfH.EndI(); TI++) {
        betweencenter << proteins.at(TI.GetKey()) << " "<< TI.GetKey() << " "<< TI.GetDat() << endl;
    }
}

//Analyzes EigenVectorCentr
void analyze_GetEigenVectorCentr(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins) {
    string output_path = get_output_path("eigenVectorCenter", organism, simple_output);
    ofstream eigvectorcenter(output_path.c_str());
    TIntFltH NIdVCfH;
    TSnap::GetEigenVectorCentr(G, NIdVCfH);

    eigvectorcenter << fixed << setprecision(PRECISION);

    for (TIntFltH::TIter TI = NIdVCfH.BegI(); TI < NIdVCfH.EndI(); TI++) {
        eigvectorcenter << proteins.at(TI.GetKey()) << " "<< TI.GetKey() << " "<< TI.GetDat() << endl;
    }
}

// Analyzes the degree inequality of the graph with Lorenz curves
void analyze_lorenz(PGraph G, int organism, bool simple_output) {
    static const int NUM_POINTS = 50;
    string input_path = get_output_path(DEGREE_DIR, organism, simple_output);
    string output_path = get_output_path(LORENZ_DIR, organism, simple_output);
    ifstream degrees(input_path.c_str(), ios_base::app);
    ofstream lorenz(output_path.c_str());

    int64 total_degree_sum = 2 * G->GetEdges();
    int64 total_nodes = G->GetNodes();
	
    lorenz << fixed << setprecision(PRECISION);

    for (int i = 1; i <= NUM_POINTS; i++) {
        int64 degree_sum = 0;
        int64 nodes = total_nodes * i / NUM_POINTS;
        string line;
        int degree, frequency;

        degrees.clear();
        degrees.seekg(0);

        while (getline(degrees, line)) {
            sscanf(line.c_str(), "%d %d", &degree, &frequency);
			
            if (nodes > frequency) {
                nodes -= frequency;
                degree_sum += frequency * degree;
            }
            else {
                degree_sum += nodes * degree;
                break;
            }
        }

        lorenz << (double)i / NUM_POINTS << " "
            << (double)degree_sum / total_degree_sum << endl;
    }
}

// Analyzes the degree centrality distribution of the graph
void analyze_degreecentrality(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins, TIntFltH &NIdDegH) {
    string output_path = get_output_path(DEGREECTR_DIR, organism, simple_output);
    ofstream degreecentrality(output_path.c_str());

    degreecentrality << fixed << setprecision(PRECISION);

    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        int NId = NI.GetId();
        double Deg = TSnap::GetDegreeCentr(G, NId);
        NIdDegH.AddDat(NId, Deg);
        degreecentrality << proteins.at(NId) <<" NId:"<<NId<< " Deg:"<< Deg << endl;
    }
}

//Analyzes the close centrality distribution of the graph
void analyze_closecentrality(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins) {
    string output_path = get_output_path("closecenter", organism, simple_output);
    ofstream closecentrality(output_path.c_str());

    closecentrality << fixed << setprecision(PRECISION);

    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        int NId = NI.GetId();
        double close = TSnap::GetClosenessCentr(G, NId);
		closecentrality << proteins.at(NId) << " "<< NId << " "<< close << endl;
    }
}
// Analyzes the normalized node betweenness centrality distribution of the graph
void analyze_betweenness(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins, TIntFltH &NIdBtwH) {
    string output_path = get_output_path(BETWEEN_DIR, organism, simple_output);
    //ofstream betweenness(output_path.c_str(), ios_base::app);
    ofstream betweenness(output_path.c_str());
    TSnap::GetBetweennessCentr(G, NIdBtwH);
    int nodes = G->GetNodes();
    long long max_nodes = ((long long)nodes - 1) * (nodes - 2) / 2;

    betweenness << fixed << setprecision(PRECISION);

    for (TIntFltH::TIter TI = NIdBtwH.BegI(); TI < NIdBtwH.EndI(); TI++) {
        betweenness << proteins.at(TI.GetKey()) << " "
            << TI.GetDat() / max_nodes << endl;
    }
}

// Analyzes the closeness centrality distribution of the graph
void analyze_closeness(PGraph G, int organism, bool simple_output,
        const vector<string> &proteins, TIntFltH &NIdCloH) {
    string output_path = get_output_path(CLOSE_DIR, organism, simple_output);
    ofstream closeness(output_path.c_str());

    closeness << fixed << setprecision(PRECISION);

    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        int NId = NI.GetId();
        double Clo = TSnap::GetClosenessCentr<PGraph>(G, NId, false);
        NIdCloH.AddDat(NId, Clo);
        closeness << proteins.at(NId) << " "
            << Clo << endl;
    }
}

// Analyzes the eigenvalues of the adjacency matrix of the graph
void analyze_adjacency(PGraph G, int organism, bool simple_output) {
    string output_path = get_output_path(ADJACENCY_DIR, organism, simple_output);
    ofstream adjacency(output_path.c_str());
    TFltV EigValV;
    TVec<TFltV> EigVecV;
    int nodes = G->GetNodes();
    TSnap::GetEigVec(G, nodes, EigValV, EigVecV);

    adjacency << fixed << setprecision(PRECISION);

    for (TFltV::TIter TI = EigValV.BegI(); TI < EigValV.EndI(); TI++) {
        adjacency << *TI << endl;
    }
}

// Analyzes fragmentation for a single node removal strategy

void node_removal(PGraph G, int organism, string strategy, bool simple_output, TIntFltH &Hash) {
    static const int NUM_POINTS = 100;
    string output_path = get_output_path(FRAGMENT_DIR, organism, simple_output);
    output_path += strategy;

    ofstream fragmentation(output_path.c_str());
    fragmentation << fixed << setprecision(PRECISION);

    int nodes = G->GetNodes();

	fragmentation<<nodes << endl;
    Hash.SortByDat(false);
	
	//Returns iterator to the first element of the hash table.
    TIntFltH::TIter TI = Hash.BegI(); 
    int num_deleted = 0;
    PGraph H = new PGraph::TObj(*G);
    int current_nodes;
	
    for (int i = 1; i <= NUM_POINTS; i++) {	
        for ( ; num_deleted < nodes * i / NUM_POINTS; TI++, num_deleted++) {  //按照比例移出
            H->DelNode(TI.GetKey());
        }
		
		//cout<< "---------enddel---------" << endl;  //add diy
        TIntPrV SccSzCnt;
		//SccSzCnt	returns a set of pairs (number of nodes in the component, number of such components)
        TSnap::GetSccSzCnt(H, SccSzCnt);
        current_nodes = H->GetNodes();
        fragmentation << (double)i / NUM_POINTS << " " << current_nodes << ", ";

        for (TIntPrV::TIter TI1 = SccSzCnt.BegI(); TI1 < SccSzCnt.EndI(); TI1++) {			
            TInt size = TI1->GetVal1();
            TInt frequency = TI1->GetVal2();
            fragmentation << size << " " << frequency << ", ";
        }
        fragmentation << endl;
    }
}

// Analyzes the fragmentation of the graph to node and edge removal
void analyze_fragmentation(PGraph G, int organism, bool simple_output,
        TIntFltH &NIdDegH, TIntFltH &NIdBtwH, TIntFltH &NIdCloH) {
    TIntFltH NIdRndH;
    TRnd random(0); 
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        NIdRndH.AddDat(NI.GetId(), random.GetUniDev());   
    }

    node_removal(G, organism, "random", simple_output, NIdRndH); // remove random
}

// Creates graph of organism's data and analyzes it
void analyze_organism(int organism, string input, bool simple_output,int rates) {
    cout << "Started organism " << organism << endl;
    vector<string> proteins;
    TIntFltH NIdDegH, NIdBtwH, NIdCloH;
	int i=0;int j=0; double sum=0; int maxcout=10000;
	PGraph H;
	
    PGraph G = create_graph(organism, input, proteins);	
	ofstream resaverage("Del1Noderes100",ios::app);  
	int edges=G->GetEdges();
	int nodes = G->GetNodes();
	int rmvnode=0;
	int firstflag=1;
	
    G = TSnap::GetMxScc(G); 
	edges=G->GetEdges();
	nodes = G->GetNodes();
	cout<<"GetMxScc nodes: "<<nodes<<", edges: "<<edges<<endl; 

	resaverage<<"====newstart maxnodes "<<nodes<<" maxedges "<<edges<<endl; 

	analyze_statistics(G, organism, simple_output); 
	analyze_kcores(G, organism, simple_output);
	analyze_degrees(G, organism, simple_output);  
	analyze_nodedegree(G, organism, simple_output,proteins);
	analyze_curvature(G, organism, simple_output, proteins);
	analyze_clustering(G, organism, simple_output, proteins);
	analyze_lorenz(G, organism, simple_output);  
	analyze_degreecentrality(G, organism, simple_output, proteins, NIdDegH);
	analyze_closecentrality(G, organism, simple_output, proteins); 
	analyze_GetBetweennessCentr(G, organism, simple_output, proteins);
	analyze_GetEigenVectorCentr(G, organism, simple_output, proteins);

	for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++)
	{
		H = new PGraph::TObj(*G); 	
		rmvnode=NI.GetId();
		H->DelNode(NI.GetId());			

		sum=0;
		for(i=0;i<maxcout;i++) //计算100次,求平均
		{
		analyze_fragmentation(H, organism, simple_output, NIdDegH, NIdBtwH, NIdCloH);
		
		double ret=res(organism, "random", false); 

		sum +=ret;
		}
		resaverage<<organism<<" rmvnode: "<<rmvnode<<" res: "<<(sum/maxcout)<<endl; 
	}

	resaverage.close();
	cout << "Finished organism " << organism << endl;
}

void NodeSortByRandom(PGraph G)
{
	TIntFltH::TIter TISortLocal;
	TRnd random(0); // 0 means seed from clock
	for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        NIdRndH11.AddDat(NI.GetId(), random.GetUniDev());   
    }
	
    NIdRndH11.SortByDat(false);
	
	TISort = NIdRndH11.BegI(); 
	
	cout<<"===sort:===="<<endl;
	TISortLocal=TISort;
	for (int num_deleted = 0; num_deleted < 44; TISortLocal++, num_deleted++) 
	{ 
		cout<<"id:"<<TISortLocal.GetKey()<<" randomseed:"<<NIdRndH11.GetDat(TISortLocal.GetKey())<<endl;		
    }
}
int beforeAnalyzeRemove(PGraph G,int removepercent)
{

	int nodes = G->GetNodes();
	int recRmvNode;
	int removenode=removepercent;   

	for (int num_deleted = 0; num_deleted < removenode; TISort++, num_deleted++) 
	{ 

		G->DelNode(TISort.GetKey());
		recRmvNode=TISort.GetKey();
		
    }
	return recRmvNode;
}
double res(int organism, string strategy,bool simple_output) {

	int Div_Count = 100;
	int NUM_POINTS;
    string input_path = get_output_path(FRAGMENT_DIR, organism, simple_output);
    input_path += strategy;
    ifstream fragmentation(input_path.c_str());//, ios_base::app
	
	string output_path = get_output_path(FRAGMENT_DIR, organism, simple_output);
    output_path += strategy;
	output_path += "hmsh";
	
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
	double ret=0;
	string x ;
	double pilog2=0;
	double pilog3=0;
	double percent=0;
	
	fragmentation.clear();
	fragmentation.seekg(0);
	
	getline(fragmentation, line);
	sscanf(line.c_str(), "%d", &NUM_POINTS);  //取节点数，逐一移出	

	while (getline(fragmentation, line)) {
		pi=0;
		pi2=0;
		pilog=0;
		pilog2=0;	
		
		index = line.find(" ",0);
		if(index < line.length())
			;
		else
			;
		x = line.substr(0,index);
		sscanf(x.c_str(), "%lf", &percent);
		
		size_t size = line.size();
	
		y= line.substr(index+1,size);

		sscanf(y.c_str(), "%d", &nodesum);
		
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

			if(NUM_POINTS-nodesum>0)
			{
				pi2=(double)(NUM_POINTS-nodesum)/NUM_POINTS;				
				pilog2 = pi2*log((double)1/NUM_POINTS);
								
			}

			pilog= pilog +pilog2;
			pilog= pilog*(-1);
			pilog= pilog/log(NUM_POINTS);   //modified shannon	
		}
		
		if((nodesum == 1)||(nodesum == 0)) //isolate node
		{

			pilog = 1;
		}		
		pilog3 +=pilog;  
	}
	
	pilog3=pilog3/Div_Count;
	ret = 1-pilog3;  

	return ret;
}
// Runs desired function in parallel
void run_in_parallel() {
    thread threads[NUM_ORG];
    ifstream org_list(ORG_LOCATION.c_str());
    string line;
    for (int i = 0; getline(org_list, line); i++) {
        if (i >= CONCURRENCY) threads[i - CONCURRENCY].join();

        int organism = stoi(line);
        threads[i] = thread(analyze_organism, organism, "", false,0);
    }
    for (int i = NUM_ORG - CONCURRENCY; i < NUM_ORG; i++) {
        threads[i].join();
    }
}

// Generates list of commands
void generate_job_list(int num_instances) {
    ifstream org_list(ORG_LOCATION.c_str());
    ofstream job_list(JOB_LIST_LOCATION.c_str());
    string line;

    while(getline(org_list, line)) {
        for (int j = 0; j <= num_instances; j++) {
            job_list << "./analyze " << line << " " << j << " "
                << EXTRACTED_DATA_DIR << "/" << line << endl;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        make_analyze_dir();
        run_in_parallel();
    }
    else if (argc > 4) {
        cout << "Use one of the following formats:" << endl
            << "[executable] to run all organisms in parallel" << endl
            << "[executable] [organism] to run [organism]" << endl
            << "[executable] [organism]" << endl
            << "[executable] [organism] [input]" << endl;
    }
    else {
        int organism = atoi(argv[1]);
        string input = (argc == 4) ? argv[2]: "";
        analyze_organism(organism, input, false,atoi(argv[3]));
    }
    return 0;
}
