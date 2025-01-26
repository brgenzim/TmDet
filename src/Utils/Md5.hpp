// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

/**
 * @file
 * @author [tGautot](https://github.com/tGautot)
 * @brief Simple C++ implementation of the [MD5 Hashing
 * Algorithm](https://en.wikipedia.org/wiki/MD5)
 * @details
 * The [MD5 Algorithm](https://en.wikipedia.org/wiki/MD5) is a
 * hashing algorithm which was designed in 1991 by [Ronal
 * Rivest](https://en.wikipedia.org/wiki/Ron_Rivest).
 *
 * MD5 is one of the most used hashing algorithm there is. Some of its
 * use cases are:
 *  1. Providing checksum for downloaded software
 *  2. Store salted password
 *
 * However MD5 has be know to be cryptographically weak for quite some
 * time, yet it is still widely used. This weakness was exploited by the
 * [Flame Malware](https://en.wikipedia.org/wiki/Flame_(malware)) in 2012
 *
 * ### Algorithm
 * First of all, all values are expected to be in [little endian]
 * (https://en.wikipedia.org/wiki/Endianness). This is especially important
 * when using part of the bytestring as an integer.
 *
 * The first step of the algorithm is to pad the message for its length to
 * be a multiple of 64 (bytes). This is done by first adding 0x80 (10000000)
 * and then only zeroes until the last 8 bytes must be filled, where then the
 * 64 bit size of the input will be added
 *
 * Once this is done, the algo breaks down this padded message
 * into 64 bytes chunks. Each chunk is used for one *round*, a round
 * breaks the chunk into 16 blocks of 4 bytes. During these rounds
 * the algorithm will update its 128 bit state (represented by 4 ints: A,B,C,D)
 * For more precisions on these operations please see the [Wikipedia
 * aritcle](https://en.wikipedia.org/wiki/MD5#Algorithm).
 * The signature given by MD5 is its 128 bit state once all rounds are done.
 * @note This is a simple implementation for a byte string but
 * some implmenetations can work on bytestream, messages of unknown length.
 */

#pragma once

#include <algorithm>  /// Used for std::copy
#include <array>      /// Used for std::array
#include <cassert>    /// Used for assert
#include <cstring>    /// Used for std::memcopy
#include <iostream>   /// Used for IO operations
#include <string>     /// Used for strings
#include <vector>     /// Used for std::vector

/**
 * @namespace hashing
 * @brief Hashing algorithms
 */
namespace hashing {
    /**
    * @struct md5
    * @brief Functions for the [MD5](https://en.wikipedia.org/wiki/MD5) algorithm
    * implementation
    */
    struct md5 {
        /**
        * @brief Rotates the bits of a 32-bit unsigned integer
        * @param n Integer to rotate
        * @param rotate How many bits for the rotation
        * @return uint32_t The rotated integer
        */
        uint32_t leftRotate32bits(uint32_t n, std::size_t rotate);

        /**
        * @brief Checks whether integers are stored as big endian or not
        * @note Taken from [this](https://stackoverflow.com/a/1001373) StackOverflow
        * post
        * @return true IF integers are detected to work as big-endian
        * @return false IF integers are detected to work as little-endian
        */
        bool isBigEndian();

        /**
        * @brief Sets 32-bit integer to little-endian if needed
        * @param n Number to set to little-endian (uint32_t)
        * @return uint32_t param n with binary representation as little-endian
        */
        uint32_t toLittleEndian32(uint32_t n);

        /**
        * @brief Sets 64-bit integer to little-endian if needed
        * @param n Number to set to little-endian (uint64_t)
        * @return uint64_t param n with binary representation as little-endian
        */
        uint64_t toLittleEndian64(uint64_t n);

        /**
        * @brief Transforms the 128-bit MD5 signature into a 32 char hex string
        * @param sig The MD5 signature (Expected 16 bytes)
        * @return std::string The hex signature
        */
        std::string sig2hex(void* sig);

        /**
        * @brief The MD5 algorithm itself, taking in a bytestring
        * @param input_bs The bytestring to hash
        * @param input_size The size (in BYTES) of the input
        * @return void* Pointer to the 128-bit signature
        */
        void* hash_bs(const void* input_bs, uint64_t input_size);

        /**
        * @brief Converts the string to bytestring and calls the main algorithm
        * @param message Plain character message to hash
        * @return void* Pointer to the MD5 signature
        */
        void* hash(const std::string& message);
    };
}  // namespace hashing

namespace Tmdet::Utils {
    struct Md5 {
        static std::string getHash(std::string source);
    };
}