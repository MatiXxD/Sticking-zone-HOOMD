#include <iostream>
#include <fstream>
#include <string>


class MakeDiagram{

    std::string s;
    std::fstream f;

    const unsigned N = 3600;
    const unsigned SPLIT_PARAM = 15;
    const unsigned FRAMES_HANDLE = 10;
    const unsigned FRAMES_SKIP = 20;

    const double LEN_X = 569.8768003947 * 2.0;

    double *x;
    double *arrN;

public:
    MakeDiagram(std::string name, int epsNum);
    ~MakeDiagram();
    void readHeader();
    void skipFrame();
    void resetN();
    void readData();
    void getN();
    void writeRes(int outNum);
};

MakeDiagram::MakeDiagram(std::string name, int epsNum){

    f.open(name.c_str());

    x = new double[N];
    arrN = new double[SPLIT_PARAM];

    for(int rep = 0 ; rep < epsNum ; rep++){
        
        // Skip frames which correspond to system's relaxation
        for(int i = 0 ; i < FRAMES_SKIP ; i++){
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
        writeRes(epsNum);

    }

}

MakeDiagram::~MakeDiagram(){

    delete[] x;
    delete[] arrN;

    f.close();

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
        if(temp < SPLIT_PARAM) arrN[temp]++;
        else arrN[SPLIT_PARAM-1]++;
    }

}

void MakeDiagram::writeRes(int outNum){
 
    std::ofstream fout;
    std::string name = "/home/ilya/MD/phaseDiagram/dump/main/dump_" + std::to_string(outNum) + ".txt";
    fout.open(name.c_str());

    for(int i = 0 ; i < SPLIT_PARAM ; i++)
        fout << i * outNum * (LEN_X / SPLIT_PARAM) << ' ' << arrN[i] / (N * SPLIT_PARAM) << '\n';

    fout.close();

}

int main(int argc, char* argv[]){

    int epsNum = 1;
    std::string name = "dump_test.txt";

    MakeDiagram(name, epsNum);

    return 0;

}
