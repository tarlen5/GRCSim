#include <highfive/highfive.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace HighFive;

bool fileExist(const std::string& filename) {
    std::ifstream file(filename.c_str());
    return file.good();
}

 // Define dataset names
std::vector<std::string> datasetNames = {"egyPrim", "egySec", "theta", "phi", "time", "thetap", "xi", "weight"};


void To_HDF5(const std::string& filename,
             const std::vector<std::vector<double>>& data) {
    try {

           // Create or open the specified group
           if (fileExist(filename))
                 {
                  File file(filename, HighFive::File::ReadWrite);
                  for (size_t i = 0; i < datasetNames.size(); ++i)  {
                  DataSet dataset = file.getDataSet(datasetNames[i]);
                  size_t currentSize = dataset.getDimensions()[1];
                  dataset.resize({1,currentSize+1});
                  dataset.select({0,currentSize},{1,1}).write(data[i]); }
                 }
           else  {
                             // Create a new file if it doesn't exist
                  File file(filename, HighFive::File::ReadWrite | HighFive::File::Create);
                   DataSpace dataspace = DataSpace({1, 1}, {100, DataSpace::UNLIMITED});
                   DataSetCreateProps props;
                   props.add(Chunking(std::vector<hsize_t>{2, 2}));
                  for (size_t i = 0; i < datasetNames.size(); ++i) {
                  DataSet dataset = file.createDataSet(datasetNames[i],dataspace,  create_datatype<double>(),props);
                  dataset.select({0, 0}, {1, 1}).write(data[i]);}
                 }
         

        // Print a message indicating the dataset and group name
        std::cout << "Data appended to " << filename  << "'" << std::endl;

    } catch (const HighFive::Exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;}
        // You may choose to throw an exception or handle the error in some other way 
}

int main() {
    // Example usage
    std::vector<std::vector<double>> Data = {
        {3.57e6},
        {5e5},
        {55.74},
        {0.57},
        {1.4e4},
        {0.6},
        {0.778},
        {0.34}
    };
    To_HDF5("histograms.h5", Data);

    return 0;
}

