#include <highfive/highfive.hpp>
#include <vector>
#include <iostream>

using namespace HighFive;

void To_HDF5(const std::string& filename,
             const std::string& groupName,
             const std::string& datasetName,
             const std::vector<std::vector<double>>& data,
             const std::string& xAxis,
             const std::string& yAxis) {
    try {
        // Create a new HDF5 file for writing
        File file(filename, HighFive::File::ReadWrite | HighFive::File::Create);

        // Create the specified group
        auto group = file.createGroup(groupName);

        // Create the specified dataset within the group
        DataSet dataset = group.createDataSet(datasetName, data);

        // Create XAxis and YAxis attributes for the dataset
        dataset.createAttribute("XAxis", std::string(xAxis));
        dataset.createAttribute("YAxis", std::string(yAxis));

        // Print a message indicating the dataset and group name
        std::cout << "Dataset '" << datasetName << "' saved to " << filename << " under group '" << groupName << "'" << std::endl;

    } catch (const HighFive::Exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        // You may choose to throw an exception or handle the error in some other way
    }
}

int main() {
    // Example usage
    std::vector<std::vector<double>> egyPrimData = {
        {1e6, 2e6, 3e6},
        {100, 150, 300}
    };

    To_HDF5("histograms.h5", "egyPrim", "Histogram", egyPrimData, "egyPrim", "Num");


    return 0;
}
