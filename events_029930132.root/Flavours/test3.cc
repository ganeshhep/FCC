#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <vector>

bool CompareColumns(const std::vector<int>& column1, const std::vector<int>& column2) {
    // Check if both columns have the same size
    if (column1.size() != column2.size()) {
        return false; // Sizes don't match, so columns cannot match
    }

    // Compare corresponding elements
    for (size_t i = 0; i < column1.size(); ++i) {
        if (column1[i] != column2[i]) {
            return false; // Found a mismatch
        }
    }

    return true; // All elements match
}

void test3()
{
    std::vector<int> col1 = {1, 2, 3, 4, 5};
    std::vector<int> col2 = {1, 2, 3, 4, 5};
    std::vector<int> col3 = {1, 2, 3, 4, 6};

    std::cout << "col1 and col2 match: " << CompareColumns(col1, col2) << std::endl;
    std::cout << "col1 and col3 match: " << CompareColumns(col1, col3) << std::endl;

}
