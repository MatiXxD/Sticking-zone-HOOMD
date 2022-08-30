#include <iostream>
#include <fstream>
#include <string>

#define EPS_NUM 13

//0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0
//1    2    3    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18

class MakeDiagram{

    std::string s;
    std::fstream f;

    const unsigned N = 3600;
    const unsigned SPLIT_PARAM = 15;
    const unsigned FRAMES_HANDLE = 10;
    const unsigned FRAMES_SKIP = 20;

    const double LEN_X = 341.9260802368 * 2.0;

    double *x;
    double *arrN;


public:
    MakeDiagram(std::string dir, std::string name);
    ~MakeDiagram();
    void readHeader();
    void skipFrame();
    void resetN();
    void readData();
    void getN();
    void writeRes(int outNum);
};

MakeDiagram::MakeDiagram(std::string dir, std::string name){

    x = new double[N];
    arrN = new double[SPLIT_PARAM];

    int skip;
    for(int rep = 0 ; rep < EPS_NUM ; rep++){
        
        // Open file to read data 
        std::string nameOut = dir + name + std::to_string(rep+1) + ".txt";
        std::cout << "Repeat: " << rep << std::endl;
        f.open(nameOut.c_str());

        // Skip frames which correspond to system's relaxation
        if(rep != 0)
            skip = FRAMES_SKIP;
        else 
            skip = FRAMES_SKIP * 10;
        for(int i = 0 ; i < skip ; i++){
            readHeader();
            skipFrame();
        }
        
        // Reset array with particles distribution
        resetN();

        // Get avg rho for current epsilon value
        for(int i = 0 ; i < FRAMES_HANDLE ; i++){
            readHeader();
            readData();
            getN();
        }
        
        // Write results in file
        writeRes(rep);

        // Close file
        f.close();

    }

}

MakeDiagram::~MakeDiagram(){

    delete[] x;
    delete[] arrN;

}

void MakeDiagram::readHeader(){

    f >> s; f >> s;
    f >> s;
    f >> s; f >> s; f >> s; f >> s;
    f >> s;
    f >> s; f >> s; f >> s; f >> s; f >> s; f >> s;
    f >> s; f >> s;
    f >> s; f >> s;
    f >> s; f >> s;
    f >> s; f >> s; f >> s; f >> s; f >> s;

}

void MakeDiagram::skipFrame(){

    for(int i = 0 ; i < N ; i++){
        f >> s; f >> s; f >> s;
    }

}

void MakeDiagram::resetN(){

    for(int i = 0 ; i < SPLIT_PARAM ; i++)
        arrN[i] = 0.0;

}

void MakeDiagram::readData(){
    
    int index;

    for(int i = 0 ; i < N ; i++){
        f >> index;
        f >> x[index];
        f >> s;
    }    

}

void MakeDiagram::getN(){

    int temp;

    for(int i = 0 ; i < N ; i++){
        temp = (x[i] + LEN_X / 2.0) / (LEN_X / SPLIT_PARAM);
        if(temp < SPLIT_PARAM) 
            arrN[temp]++;
        else 
            arrN[SPLIT_PARAM-1]++;
    }

}

void MakeDiagram::writeRes(int outNum){
 
    std::ofstream fout;
    std::string name = "/home/ilya/MD/phaseDiagram/dump/main2/dump_" + std::to_string(outNum+1) + ".txt";
    fout.open(name.c_str());

    for(int i = 0 ; i < SPLIT_PARAM ; i++)
        fout << i * (LEN_X / SPLIT_PARAM) << ' ' << arrN[i] / (N * FRAMES_HANDLE) << '\n';

    fout.close();

}

int main(int argc, char* argv[]){

    std::string dir = "/home/ilya/MD/phaseDiagram/dump/hoomdDump2/";
    std::string name = "d_v";

    MakeDiagram(dir, name);

    return 0;

}
