/*****************************************************************
 *
 *  Code: Connectivity Driven Early Floorplaning
 *  Author: A0230289
 *
 *****************************************************************
 *
 *  Usage: 
 *  ******
 *
 *  $0 PARTITION_AREA_FILE PARTITION_CONNECTIVITY_FILE OUTPUT_FILE
 *
 *  Format of PARTITION_AREA_FILE:
 *  ******************************
 *
 *  NUMBER_OF_PARTITIONS
 *  for each of these partitions:
 *      PARTITION_NAME SOFT_AREA HARD_MACRO_COUNT
 *      for each of these hard macros:
 *          HARDMACRO_WIDTH HARD_MACRO_HEIGHT
 *  
 *  Eg:
 *  3
 *  BSS 122241 2
 *  345 345
 *  120 100
 *  MSS 1234 0
 *  TOP 421231 0
 *
 *  Format of PARTITION_CONNECTIVITY_FILE:
 *  **************************************
 *  
 *  NUMBER_OF_CONNECTIONS
 *  for each connection:
 *      FROM_INDEX TO_INDEX CONNECTION_COUNT
 *
 *  Eg:
 *  3
 *  BSS TOP 200
 *  MSS BSS 500
 *  TOP MSS 420
 *
 */


#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <stack>
#include <queue>
#include "./bin_pack/MaxRectsBinPack.h"
#include "./bin_pack/MaxRectsBinPack.cpp"
#include "./bin_pack/Rect.cpp"

using namespace std;

//Constants
double MIN_MICRON_PRECISION = 1;
double MAX_MICRON_PRECISION = 50;
double MINIMUM_PARTITION_EDGE_SIZE = 50;
double AREA_COMPENSATION_FACTOR = 1e6;
double WIRE_SPACING = 0.15;

//Meta Parameters
int GENERATION_COUNT = 10;
int POOL_SIZE = 30;
double FITNESS_FACTOR = 0.2;
double MUTATION_FACTOR = 0.1;
double CROSSOVER_FACTOR = 0.3;
double DIVERSITY_FACTOR = 0.1;

//Structs in 2D Bin Packing
struct Rectangle {
    double width, height;
    //We want to sort in decreasing. MAXRECTS-BSSF-DESCSS.
    bool operator<(Rectangle other) const {
        return  min(width,height) > min(other.width,other.height);
    }
};
Rectangle CORE_BOX;

struct Coordinate {
    double x, y;
};

//Partition Class
class Partition {
    string name;
    double total_area;
    double max_width;
    double max_height;
    double soft_area;
    vector<Rectangle> hard_macros;
    public:
    Partition(string name, double soft_area, vector<Rectangle> hard_macros) {
        assert(soft_area>=0);
        this->total_area = -1;
        this->max_width = -1;
        this->max_height = -1;
        this->name = name;
        this->soft_area = soft_area;
        for (int i=0; i<hard_macros.size(); i++) {
            assert(hard_macros[i].width > 0);
            assert(hard_macros[i].height > 0);
        }
        this->hard_macros = hard_macros;
    }
    double totalArea() {
        if (this->total_area == -1) {
            this->total_area = this->soft_area;
            for (int i=0; i<this->hard_macros.size();i++) {
                Rectangle hm = this->hard_macros[i];
                this->total_area += hm.width * hm.height;
            }
        }
        return this->total_area;
    }
    double maxWidth() {
        if (this->max_width == -1 && this->hard_macros.size() > 0) {
            for (int i=0; i<this->hard_macros.size(); i++) {
                if (this->max_width == -1 || this->hard_macros[i].width > this->max_width) {
                    this->max_width = this->hard_macros[i].width;
                }
            }
        }
        return max(this->max_width,MINIMUM_PARTITION_EDGE_SIZE);
    }
    double maxHeight() {
        if (this->max_height == -1 && this->hard_macros.size() > 0) {
            for (int i=0; i<this->hard_macros.size(); i++) {
                if (this->max_height == -1 || this->hard_macros[i].height > this->max_height) {
                    this->max_height = this->hard_macros[i].height;
                }
            }
        }
        return max(this->max_height,MINIMUM_PARTITION_EDGE_SIZE);
    }
    vector<Rectangle> getHardMacros() {
        return this->hard_macros;
    }
    double getSoftArea() {
        return this->soft_area;
    }
    double getHardArea() {
        return this->totalArea() - this->soft_area;
    }

};
vector<Partition> partitions;
vector<pair<string,int> > partition_name_to_index;

//Connection Struct
struct Connection {
    int from, to, count;
};
vector<Connection> connections;

void readInput(string partition_area_file, string partition_connection_file) {
    cout << "Reading Inputs." << endl;
    ifstream input;
    input.open(partition_area_file.c_str());
 
    //Line 1:
    int partition_count;
    input >> partition_count;

    assert(partition_count > 2);
    partitions = vector<Partition>();
    partition_name_to_index = vector<pair<string,int> >();
    for (int i=0; i<partition_count; i++) {
        //Line 2*: NAME SOFT_AREA HARD_MACRO_COUNT
        string partition_name;
        double soft_area;
        int hard_macro_count;
        input >> partition_name >> soft_area >> hard_macro_count;
        vector<Rectangle> hard_macros = vector<Rectangle>();
        for (int j=0; j<hard_macro_count;j++) {
            Rectangle hard_macro;
            input >> hard_macro.width >> hard_macro.height;
            hard_macros.push_back(hard_macro);
        }
        partitions.push_back(Partition(partition_name,soft_area,hard_macros));
        partition_name_to_index.push_back(make_pair(partition_name,i));
    }   
    input.close();
    cout << "Read " << partitions.size() << " partitions." << endl; 

    input.open(partition_connection_file.c_str());
    int connection_count;
    input >> connection_count;
    connections = vector<Connection>();
    for (int i=0; i<connection_count; i++) {
        Connection connection;
        string from, to;
        input >> from >> to >> connection.count;
        int from_index=-1;
        for (int k=0; k<partition_name_to_index.size(); k++) {
            if (partition_name_to_index[k].first == from) {
                from_index=k;
                break;
            }
        }
        assert(from_index>=0);
        int to_index=-1;
        for (int k=0; k<partition_name_to_index.size(); k++) {
            if (partition_name_to_index[k].first == to) {
                to_index=k;
                break;
            }
        }
        assert(to_index>=0);
        assert(from_index != to_index);
        connection.from = from_index;
        connection.to   = to_index;
        connections.push_back(connection);
    }
    input.close();
    double total_area = 0;
    for (int i=0; i<partitions.size(); i++) {
        total_area += partitions[i].totalArea();
    }
    CORE_BOX.width = CORE_BOX.height = sqrt(total_area);
}

/*
 *  Interface to MAXRECTS-BSSF-DESCSS
 *  *********************************
 *  
 *  I have made a few changes to the original code available on GitHub:
 *  https://github.com/juj/RectangleBinPack
 *
 *  The original code was implemented with the 'int' data type.
 *  Changed this->to 'double'.
 *
 *  Important:
 *  **********
 *
 *  Only made changes to 'BSSF' code so if you want to expriment with other
 *  heuristics, you will have to modify the 'MaxRectsBinPack.cpp' file in
 *  the same directory and re-compile. Also: no 90 deg fliping by default.
 *
 */

bool can2DPack (Rectangle container, vector<Rectangle> contents, bool verbose=false) {
    //DESCSS
    sort(contents.begin(), contents.end());
    //Namespace rbp contains the functions we need.
    rbp::MaxRectsBinPack bin;
    bin.Init(container.width, container.height,false);
    for (int i=0; i<contents.size(); i++) {
        //BSSF is default insertion scheme.
        rbp::Rect packedRect = bin.Insert(contents[i].width, contents[i].height, rbp::MaxRectsBinPack::RectBestShortSideFit);
        if (packedRect.height == 0) {
            return false;
        }
        if (verbose) {
            double x1 = packedRect.x;
            double y1 = packedRect.y;
            double x2 = x1 + packedRect.width;
            double y2 = y1 + packedRect.height;
            cout << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
        }
    }
    return true;
}

//Lookup Table
vector<vector<bool> > legalPartitionWidthsLUT;
vector<double> legalPartitionWidthLUTPrecisions;

void initializeLookupTables() {
    //First set precision of evaluation i.e. the delta-w.
    //This depends on how large the area is.
    //We currently don't carea about the ratio of hard and soft area.
    cout << "Initializing Lookup Tables." << endl;
    legalPartitionWidthLUTPrecisions = vector<double>(partitions.size());
    for (int i=0; i<partitions.size(); i++) {
        double width_precision = max(MIN_MICRON_PRECISION,min(partitions[i].totalArea()/AREA_COMPENSATION_FACTOR,MAX_MICRON_PRECISION));
        legalPartitionWidthLUTPrecisions[i] = width_precision;
    }
    //Populate LUT.
    legalPartitionWidthsLUT = vector<vector<bool> >(partitions.size());
    for (int i=0; i<partitions.size(); i++) {
        cout << "Initializing LUT for Partition #" << (i+1) << "." << endl;
        legalPartitionWidthsLUT[i] = vector<bool>();

        int count = 0;
        int number_of_entries =  (partitions[i].totalArea() / partitions[i].maxHeight() - partitions[i].maxWidth()) / legalPartitionWidthLUTPrecisions[i];

        for (double width = partitions[i].maxWidth(); width < partitions[i].totalArea() / partitions[i].maxHeight(); width += legalPartitionWidthLUTPrecisions[i]) {
            double implied_height = partitions[i].totalArea() / width;
            Rectangle container;
            container.width = width, container.height = implied_height;
            legalPartitionWidthsLUT[i].push_back(can2DPack(container,partitions[i].getHardMacros()));
            ++count;
            cout << count << "/" << number_of_entries << " widths computed.\r";
        }
        cout << endl;
    }
    cout << "Lookup Tables Initialized." << endl;
}

bool isPartitionWidthLegal(int partition_index, double width) {
    int i = partition_index;
    //These LUT values dont exist as they don't need to be computed.
    if (width > partitions[i].totalArea() / partitions[i].maxHeight()) return false;
    //Widths before the hard macro maxWidths were also not computed.
    double effective_width = width - partitions[i].maxWidth();
    if (effective_width < 0) return false;
    //Round down to nearest computed legal value.
    int LUT_INDEX = effective_width / legalPartitionWidthLUTPrecisions[i];
    return legalPartitionWidthsLUT[i][LUT_INDEX];
}

//Partition + Aspect Ratio + Bottom Left Coordinate in Floorplan = Partition Box
class PartitionBox {
    Coordinate bottom_left_coordinate;
    Rectangle bounding_box;
    int partition_index;
    public:
    //Specify valid aspect ratios!
    PartitionBox(int partition_index, Rectangle bounding_box, Coordinate bottom_left_coordinate) {
        this->partition_index = partition_index;
        this->bounding_box = bounding_box;
        this->bottom_left_coordinate = bottom_left_coordinate;
    }
    PartitionBox() {
        this->partition_index=-1;   
    }
    bool hasLegalDimensions() {
        return isPartitionWidthLegal(this->partition_index, this->bounding_box.width);
    }
    int getPartitionIndex() {
        return this->partition_index;
    }
    Coordinate getOrigin() {
        return this->bottom_left_coordinate;
    }
    Rectangle getBox() {
        return this->bounding_box;
    }
};

//Abstraction to GA
class Gene {
    public:
    double fitness;
    vector<int> NPE;
    vector<PartitionBox> DNA;
    Gene(vector<int> NPE, vector<PartitionBox> DNA) {
        this->NPE = NPE;
        this->DNA = DNA;
    }
    bool operator == (const Gene& gene) const {
        if (this->NPE != gene.NPE) {
            return false;
        }
        return true;
    }
};

double getConnectedLength(PartitionBox partition_box1, PartitionBox partition_box2) {
    //Atmost one connected edge with current architecture.
    //PartitionBox1
    double p1_x1 = partition_box1.getOrigin().x;
    double p1_y1 = partition_box1.getOrigin().y;
    double p1_x2 = p1_x1 + partition_box1.getBox().width;
    double p1_y2 = p1_y1 + partition_box1.getBox().height;
    //PartitionBox2
    double p2_x1 = partition_box2.getOrigin().x;
    double p2_y1 = partition_box2.getOrigin().y;
    double p2_x2 = p2_x1 + partition_box2.getBox().width;
    double p2_y2 = p2_y1 + partition_box2.getBox().height;
    //If not touching
    if (p1_x1 > p2_x2 || p2_x1 > p1_x2 || p1_y1 > p2_y2 || p2_y1 > p1_y2) return 0;
    double o_x1 = max(p1_x1,p2_x1);
    double o_y1 = max(p1_y1,p2_y1);
    double o_x2 = min(p1_x2,p2_x2);
    double o_y2 = min(p1_y2,p2_y2);
    double delta_x = o_x2 - o_x1;
    double delta_y = o_y2 - o_y1;
    //Account for rounding errors.
    assert(!(delta_x > 1e-8 && delta_y > 1e-8));
    return delta_x + delta_y;
}

//TODO:
//Prepare this for multithreading
double evaluateGeneFitness(Gene gene) {
    //Given an NPE and the partition area vector partitions[i].totalArea().
    //The dimensions of the partitions are fixed. Find these dimensions.
    //Primary Requirement.
    for (int i=0; i < gene.DNA.size(); i++) {
        //If the dimensions are illegal for any partition box, the floorplan is illegal.
        if (!gene.DNA[i].hasLegalDimensions())  return 0;
    }
    //Compute fitness based on connectivity.
    double fitness = 0;
    double required_connected_length = 0;
    double actual_connected_length = 0;    
    for (int i=0; i < connections.size(); i++) {
        Connection connection = connections[i];
        required_connected_length += 1.0 * connection.count * WIRE_SPACING;
        //Effective actual connection length is at most the required connected length.
        actual_connected_length += min(1.0 * connection.count * WIRE_SPACING, getConnectedLength(gene.DNA[connection.from],gene.DNA[connection.to]) * WIRE_SPACING);
    }
    fitness = actual_connected_length / required_connected_length;
    return fitness;
}

//Not Exactly Uniform Random generator between 0,1.
//Unfortunately C++11 is not supported on RedHat 4.4.7
double uniform() {
    return ((double) rand() / (RAND_MAX));
}

int uniform(int m) {
    return (m * uniform());
}

vector<int> randomPermutation(int m) {
    vector<int> permutation = vector<int>();
    for (int i=0;i<m;i++) {
        permutation.push_back(i);
    }
    //Knuth shuffle
    for (int i=0; i<m-1; i++) {
        int j = i + uniform(m-i);
        int temp = permutation[i];
        permutation[i] = permutation[j];
        permutation[j] = temp;
    }
    return permutation;
}

vector<int> randomNPE() {
    int size = partitions.size();
    vector<int> number_component = randomPermutation(size);
    int number_index = 0;
    vector<int> NPE = vector<int>();
    int number_of_operands = 0;
    for (int i=0; i<2*size-1; i++) {
        //First two always numbers.
        if (i < 2 || number_of_operands < 2 || (number_index < size && uniform() < 0.5)) {
            NPE.push_back(number_component[number_index]);
            number_index++;
            number_of_operands++;
        } else if (NPE[i-1] == -1) {
            NPE.push_back(-2);
            number_of_operands--;
        } else if (NPE[i-1] == -2) {
            NPE.push_back(-1);
            number_of_operands--;
        } else if (uniform() < 0.5) {
            NPE.push_back(-1);
            number_of_operands--;
        } else {
            NPE.push_back(-2);
            number_of_operands--;
        }
    }
    assert(number_of_operands == 1);
    return NPE;
}

struct NPE_NODE {
    NPE_NODE * parent;
    int operation_or_index;
    double area;
    Rectangle bounding_box;
    Coordinate bottom_left_coordinate;
    NPE_NODE * firstOperand;
    NPE_NODE * secondOperand;
};

NPE_NODE * NPEvectorToTree(vector<int> NPE) {
    stack <NPE_NODE *> operand_stack;
    for (int i=0; i<NPE.size(); i++) {
        NPE_NODE * item = new NPE_NODE();
        if (NPE[i] >=0) {
            //If a partition.
            item->firstOperand = 0;
            item->secondOperand = 0;
            item->operation_or_index = NPE[i];
            item->parent = 0;
            item->area = partitions[NPE[i]].totalArea();
            operand_stack.push(item);
        } else {
            //Operator.
            item->secondOperand = operand_stack.top();
            operand_stack.pop();
            item->firstOperand = operand_stack.top();
            operand_stack.pop();
            item->firstOperand->parent = item;
            item->secondOperand->parent = item;
            item->operation_or_index = NPE[i];
            item->area = -1;
            item->parent = 0;
            operand_stack.push(item);
        }
    }
    //Only the top operand should be left at the end.
    assert(operand_stack.size() == 1);
    return operand_stack.top();
}

//Assign Area To All Operators
NPE_NODE * assignAreas(NPE_NODE * root) {
    //If operator
    if (root->area == -1) {
        assignAreas(root->firstOperand);
        assignAreas(root->secondOperand);
        root->area = root->firstOperand->area + root->secondOperand->area;
    }
    return root;
}

//Assign Dimensions To All Partitions
NPE_NODE * assignDimensions(NPE_NODE * root, Rectangle bounding_box, Coordinate bottom_left_coordinate) {
    root->bounding_box = bounding_box;
    root->bottom_left_coordinate = bottom_left_coordinate;
    //-2 => Top Bottom //-1 => Left Right
    if (root->operation_or_index == -2) {
        Rectangle top_half = bounding_box;
        Rectangle bottom_half = bounding_box;

        double top_area = root->firstOperand->area;
        double total_area = root->area;

        top_half.height *= top_area/total_area;
        bottom_half.height -= top_half.height;

        Coordinate top_half_coordinate = bottom_left_coordinate;
        top_half_coordinate.y += bottom_half.height;
        assignDimensions(root->firstOperand,top_half,top_half_coordinate);
        assignDimensions(root->secondOperand,bottom_half,bottom_left_coordinate);
    } else if (root->operation_or_index == -1) {
        Rectangle left_half = bounding_box;
        Rectangle right_half = bounding_box;

        double left_area = root->firstOperand->area;
        double total_area = root->area;

        left_half.width *= left_area/total_area;
        right_half.width -= left_half.width;

        Coordinate right_half_coordinate = bottom_left_coordinate;
        right_half_coordinate.x += left_half.width;
        assignDimensions(root->firstOperand,left_half,bottom_left_coordinate);
        assignDimensions(root->secondOperand,right_half,right_half_coordinate);
    }
    return root;
}

//NPE Tree to Partition Box
vector<PartitionBox> NPETreeToPartitonBoxVector (NPE_NODE * root) {
    vector<PartitionBox> partition_boxes = vector<PartitionBox>(partitions.size());
    queue<NPE_NODE*> q;
    q.push(root);
    while (!q.empty()) {
        NPE_NODE * front = q.front();
        q.pop();
        if (front->operation_or_index < 0) {
            //Operand
            q.push(front->firstOperand);
            q.push(front->secondOperand);
        } else {
            //Partition
            int partition_id = front->operation_or_index;
            Rectangle bounding_box = front->bounding_box;
            Coordinate bottom_left_coordinate = front->bottom_left_coordinate;
            partition_boxes[partition_id] = PartitionBox(partition_id, bounding_box, bottom_left_coordinate);
        }
    }
    return partition_boxes;
}

void printNPE(vector<int> NPE) {
    for (int i=0; i<NPE.size(); i++) {
        cout << NPE[i];
    }
    cout << endl;
}

vector<PartitionBox> NPEtoPartitionBoxes(vector<int> NPE) {
    NPE_NODE * NPE_ROOT = NPEvectorToTree(NPE);
    assignAreas(NPE_ROOT);
    Coordinate origin;
    origin.x = 0;
    origin.y = 0;
    assignDimensions(NPE_ROOT, CORE_BOX, origin);
    vector <PartitionBox> sizedPartitions = NPETreeToPartitonBoxVector(NPE_ROOT);
    return sizedPartitions;
}

//A random NPE + Partitions -> Partition Boxes -> Floorplan
Gene createGene() {
    //Randomly create an NPE
    vector<int> NPE = randomNPE();
    Gene gene = Gene(NPE, NPEtoPartitionBoxes(NPE));
    return gene;
}

//Mutation may create an invalid NPE by bringing an operator near the start.
bool invalidNPE(vector<int> NPE) {
    int operandCount = 0;
    for (int j=0; j<NPE.size(); j++) {
        if (NPE[j] < 0) {
            operandCount--;
        } else {
            operandCount++;
        }
        if (operandCount == 0) return true;
    }
    return false;
}

Gene mutateGene(Gene gene, double mutation_probability) {
    vector<int> backupNPE = gene.NPE;
    //Single structural change
    if (uniform() < mutation_probability) {
        bool valid = false;
        while (!valid) {
            int i=0, j=0;
            while (i == j) {
                int rand1 = uniform(gene.NPE.size());
                int rand2 = uniform(gene.NPE.size());
                i = min(rand1,rand2);
                j = max(rand1,rand2);
            }
            //Make sure first two indices are partitions and last one is operator.
            if ((i == 0 || i == 1) && gene.NPE[j] < 0) continue;
            if (j == gene.NPE.size() - 1 && gene.NPE[i] >=0) continue;       
    
            //Swap two indices.
            int temp = gene.NPE[i];
            gene.NPE[i] = gene.NPE[j];
            gene.NPE[j] = temp;
            //Convert PE to NPE.
            for (int k=1; k<gene.NPE.size(); k++) {
                if (gene.NPE[k] < 0 && gene.NPE[k] == gene.NPE[k-1]) {
                    if (gene.NPE[k-1] == -1) {
                        gene.NPE[k] = -2;
                    } else {
                        gene.NPE[k] = -1;
                    }   
                }
            }
            if (invalidNPE(gene.NPE)) {
                gene.NPE = backupNPE;
                continue;
            }
            valid = true;
        }
        //Re compute the area breakup.
        gene.DNA = NPEtoPartitionBoxes(gene.NPE);
    }
    return gene;
}

vector<int> NPEStructure(vector<int> NPE) {
    vector<int> structure = vector<int>();
    for (int i=0; i<NPE.size(); i++) {
        if (NPE[i] >= 0) {
            structure.push_back(-3);
        } else {
            structure.push_back(NPE[i]);
        }
    }   
    return structure;
}

//The breed function
Gene breedGenes(Gene parentA, Gene parentB) {
    vector<int> inChild = vector<int>(partitions.size());
    //Take the structure of parentA.
    vector<int> childNPE = NPEStructure(parentA.NPE);
    //A random slice of parentA.
    int i=0, j=0;
    while (i==j) {
        int rand1 = uniform(parentA.NPE.size());
        int rand2 = uniform(parentA.NPE.size());
        i = min(rand1,rand2);
        j = max(rand1,rand2);
    }
    for (int k=i; k<=j; k++) {
        childNPE[k] = parentA.NPE[k];
        if (parentA.NPE[k] >= 0) {
            inChild[parentA.NPE[k]] = 1;
        }
    }
    
    //Insert from j+1 to i-1 in child NPE.
    //Iterate in Parent B's j+1 to j.
    bool started = false;
    int child_index = (j+1)%parentA.NPE.size();
    for (int k=(j+1)%parentA.NPE.size(); !started || k != (j+1)%parentA.NPE.size(); k=(k+1)%parentA.NPE.size()){
        started = true;
        //Skip operators and partitions already present in children.
        if (parentB.NPE[k] < 0 || inChild[parentB.NPE[k]] ) continue;
        //Find next empty location.
        while (childNPE[child_index] != -3) {
            child_index = (child_index + 1) % parentA.NPE.size();
        }
        childNPE[child_index] = parentB.NPE[k];
        inChild[parentB.NPE[k]] = 1;
    }

    //Don't want things to break.
    for (int i=0; i<childNPE.size(); i++) {
        assert(childNPE[i]>=-2);
    }
    return Gene(childNPE, NPEtoPartitionBoxes(childNPE));
}

//Single Thread
//Only use for the best gene of the generation.
void visualizeGene(Gene gene) {
    for (int i=0; i<gene.DNA.size(); i++) {
        cout << "Partiiton " << (i+1) << ":" << endl;
        cout << "Coordinates: " << endl;
        //Partition Dimensions
        double x1 = gene.DNA[i].getOrigin().x;
        double y1 = gene.DNA[i].getOrigin().y;
        double x2 = x1 + gene.DNA[i].getBox().width;
        double y2 = y1 + gene.DNA[i].getBox().height;
        cout << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
        cout << "Area: " << endl;
        cout << "Soft  " << partitions[gene.DNA[i].getPartitionIndex()].getSoftArea() << endl;
        cout << "Hard  " << partitions[gene.DNA[i].getPartitionIndex()].getHardArea() << endl;
        cout << "Total " << partitions[gene.DNA[i].getPartitionIndex()].totalArea() << endl;
        //Validate Approximation Earlier To Actual Width and Height
        cout << "Hard Macro Coodinates (Relative):" << endl;
        vector<Rectangle> hard_macros = partitions[gene.DNA[i].getPartitionIndex()].getHardMacros();
        if (hard_macros.size() > 0 && !can2DPack(gene.DNA[i].getBox(), hard_macros, true)) {
            cout << "ERROR: Approximation Invalid! Please reduce MIN_MICRON_PRECISION, AREA_COMPENSATION_FACTOR and re-compile." << endl;
        }
        cout << endl;   
    }
}

class GenePool {
    vector<Gene> genes;
    int pool_size;
    vector <double> gene_fitness_vector;
    bool initialized;
    //Fitness, Index
    vector<pair<double,int> > decreasing_fitness_vector;
    void evaluateFitness() {
        this->gene_fitness_vector = vector<double>(genes.size());
        this->decreasing_fitness_vector = vector<pair<double,int> >(genes.size());
        //Room for multithreading
        for (int i=0; i<genes.size();i++) {
            this->gene_fitness_vector[i] = evaluateGeneFitness(genes[i]);
            this->decreasing_fitness_vector[i] = make_pair(this->gene_fitness_vector[i],i);
        }
        //Single threaded sort
        sort(this->decreasing_fitness_vector.begin(), this->decreasing_fitness_vector.end(), greater<pair<double,int> >());
    }
    public:
    void initialize() {
        for (int i=0; i<this->pool_size; i++) {
            this->genes.push_back(createGene());
        }
        this->evaluateFitness();
        this->initialized = true;
    }
    GenePool(int size) {
        this->pool_size = size;
        this->genes = vector<Gene>();
        this->initialized = false;
    }
    Gene selectGene() {
        assert(this->initialized);
        return genes[decreasing_fitness_vector[uniform(FITNESS_FACTOR * genes.size())].second];
    }
    GenePool nextGeneration() {
        assert(this->initialized);
        assert(this->genes.size() > 0);
        //Breed the genes in this generation.
        //Consider meta parameters.
        GenePool gene_pool = GenePool(this->pool_size);
        //1. Retain the fittest.
        gene_pool.genes.push_back(genes[decreasing_fitness_vector[0].second]);
        //2. Crossover factor * this->gene_pool.size()  of the fittest genes make it to the next generation.
        for (int j = 0; j < CROSSOVER_FACTOR * this->genes.size(); j++){
            gene_pool.genes.push_back(this->selectGene());
        }
        //3. Diversity factor  * this->gene_pool.size() of new Genes created
        for (int j = 0; j < DIVERSITY_FACTOR * this->genes.size(); j++){
            gene_pool.genes.push_back(createGene());
        }
        //4. Fill remainder of gene pool with children of current generation.
        while (gene_pool.genes.size() < this->genes.size()) {
            Gene parentA = this->selectGene();
            Gene parentB = this->selectGene();
            double effective_mutation_probability = MUTATION_FACTOR;
            if (parentA == parentB) {
                effective_mutation_probability = 1;
            }
            Gene childA = mutateGene(breedGenes(parentA, parentB), effective_mutation_probability);
            Gene childB = mutateGene(breedGenes(parentB, parentA), effective_mutation_probability);
            gene_pool.genes.push_back(childA);
            if (gene_pool.genes.size() < this->genes.size()) gene_pool.genes.push_back(childB);
        }
        gene_pool.initialized = true;
        gene_pool.evaluateFitness();
        return gene_pool;
    }
    double generationFitness() {
        assert(this->initialized);
        return decreasing_fitness_vector[0].first;
    }
    void visualizeFittestGene() {
        assert(this->initialized);
        visualizeGene(genes[decreasing_fitness_vector[0].second]);
    }
};

int main(int argc, char ** argv) {
    //Read Inputs
    cout << "EFP Init" << endl;
    string partition_area_file = argv[1];
    string partition_connection_file = argv[2];
    readInput(partition_area_file, partition_connection_file);
    initializeLookupTables();
    //Genetic Algorithm
    cout << "Starting Genetic Algorithm" << endl;
    GenePool population = GenePool(POOL_SIZE);
    int generation = 0;
    population.initialize();
    while (generation < GENERATION_COUNT) {
        cout << "Generation: " << generation;
        cout << " Peak Fitness: " << population.generationFitness() << endl;
        //Validate Approximation and Visualize
        population = population.nextGeneration();
        generation++;
    }
    //Output is the last visualized Gene.
    population.visualizeFittestGene();
    return 0;
}
