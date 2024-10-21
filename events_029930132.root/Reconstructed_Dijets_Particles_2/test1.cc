#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <vector>
#include <string>

std::string GetFirstElement(const std::string &entry) {
    auto start = entry.find_first_of('[') + 1;
    auto end = entry.find_first_of(' ', start);

if (end == std::string::npos) {
        end = entry.find_first_of(']', start);
    }

    return entry.substr(start, end - start);
}

void test1()
{

std::vector<std::string> column1 = {
        "[d d ubar ubar ubar ubar ubar sbar sbar sbar sbar sbar]",
        "[d ubar]",
        "[d u s s s]",
        "[dbar u u u u u u u]",
        "[g]",
        "[g]",
        "[d ubar]"
    };

for (size_t i = 0; i < column1.size(); i++) {
        std::string targetElement = GetFirstElement(column1[i]);
std::cout << targetElement << std::endl;
 }

}
