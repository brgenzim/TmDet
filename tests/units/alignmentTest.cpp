#include <iostream>
#include <string>
#include <filesystem>
#include <Utils/Alignment.hpp>

void assertTrue(std::string testDescription, bool condition, int lineNumber);
std::string fileName;

int main() {
    fileName = std::filesystem::path(__FILE__).filename();

    // Test case 1
    {
        vector<string> query = { "MET", "ALA", "ARG", "GLY" };
        vector<string> target = { "ILE", "MET", "ALA", "ARG", "GLY" };
        Tmdet::Utils::Alignment::alignSequences(query, target);
    }

    return 0;
}

void assertTrue(std::string testDescription, bool condition, int lineNumber) {
    std::cout << (condition ? "Passed: " : "Failed: ") << testDescription;
    if (!condition) {
        std::cout << " (at line " << fileName << ":" << lineNumber << ")";
    }
    std::cout << std::endl;
}
