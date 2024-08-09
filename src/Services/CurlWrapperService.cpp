#include <iostream>
#include <curl/curl.h>
#include <Services/CurlWrapperService.hpp>

/**
 * Minimalistic wrapper for curl to perform
 * API call and get non-binary data.
 */
namespace Tmdet::Services::CurlWrapperService {

    static size_t writeBufferCallback(void* contents, size_t size, size_t chunkSize, std::string* bufferPtr);

    std::string apiCall(std::string url, Status& resultCode) {
        auto curl = curl_easy_init();
        if (curl == NULL) {
            std::cerr << "curl init failed" << std::endl;
            return "";
        }
        std::string responseBodyBuffer;
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeBufferCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &responseBodyBuffer);
        curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
        auto curlCode = curl_easy_perform(curl);

        // Some error handling
        if (curlCode != CURLE_OK) {
            resultCode = Status::Error;
            std::cerr << "HTTP Request to '" << url << "' failed: " <<
                curl_easy_strerror(curlCode) << std::endl;
        } else {
            resultCode = Status::Ok;
        }
        curl_easy_cleanup(curl);

        return responseBodyBuffer;
    }

    size_t writeBufferCallback(void* contentChunk, size_t size, size_t chunkSize, std::string* bufferPtr) {
        auto currentSize = size * chunkSize;
        bufferPtr->append(static_cast<char*>(contentChunk), currentSize);
        return currentSize;
    }
}
