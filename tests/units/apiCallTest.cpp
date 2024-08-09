#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>

#include <Services/CurlWrapperService.hpp>

void assertTrue(std::string testDescription, bool condition, int lineNumber);

std::string fileName;

int main() {

    fileName = std::filesystem::path(__FILE__).filename();
    std::string testDescription;
    Tmdet::Services::CurlWrapperService::Status status;

    // Test case 1
    {
        // Fetch data
        std::string url = "https://httpbin.org/anything?fruit=apple";
        std::string result = Tmdet::Services::CurlWrapperService::apiCall(url, status);
        auto isSuccess = result.find(R"("fruit": "apple")");
        testDescription = "request fetches JSON data";
        assertTrue(testDescription, isSuccess, __LINE__);
    }

    // Test case 2
    {
        // Handle HTTP status 404
        std::string url = "https://unitmp.org/error-with-404";
        Tmdet::Services::CurlWrapperService::apiCall(url, status);
        testDescription = "response returns HTTP error status";
        auto isFailed = (status == Tmdet::Services::CurlWrapperService::Status::Error);
        assertTrue(testDescription, isFailed, __LINE__);
    }

    // Test case 3
    {
        // Handle HTTP status 500
        std::string url = "https://httpbin.org/status/500";
        Tmdet::Services::CurlWrapperService::apiCall(url, status);
        testDescription = "response returns HTTP error status";
        auto isFailed = (status == Tmdet::Services::CurlWrapperService::Status::Error);
        assertTrue(testDescription, isFailed, __LINE__);
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
