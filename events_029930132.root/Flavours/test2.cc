#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <vector>
#include <string>

bool ContainsSandSbar(const std::string &entry) {
    return (entry.find(" s ") != std::string::npos) || (entry.find(" sbar ") != std::string::npos) ||
           (entry.find("sbar]") != std::string::npos) || (entry.find("[s ") != std::string::npos) || (entry.find("s]") != std::string::npos);
}

void test2()
{

std::vector<std::string> column1 = {
        "[d d ubar ubar ubar ubar ubar sbar sbar sbar sbar sbar]",
        "[d ubar]",
        "[d u sbar]",
        "[dbar u u u u u u u]",
        "[g s]",
        "[g]",
        "[d s ubar]"
    };

for (size_t i = 0; i < column1.size(); i++) {
        bool targetElement = ContainsSandSbar(column1[i]);
std::cout << targetElement << std::endl;
 }

int totalCount = 0; // Initialize count

    // Loop over each entry in the column
    for (const auto &entry : column1) {
        // Check if "s" or "sbar" is present
        if (ContainsSandSbar(entry)) {
            totalCount += 1; // Add 1 if "s" or "sbar" is found
        }
    }

    std::cout << "Total occurrences of 's' or 'sbar': " << totalCount << std::endl;

}
