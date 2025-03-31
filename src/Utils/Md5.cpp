#include <string>
#include <cstdint>
#include <Utils/Md5.hpp>

namespace hashing {
    uint32_t md5::leftRotate32bits(uint32_t n, std::size_t rotate) {
        return (n << rotate) | (n >> (32 - rotate));
    }

    bool md5::isBigEndian() {
        union {
            uint32_t i;
            std::array<char, 4> c;
        } bint = {0x01020304};

        return bint.c[0] == 1;
    }

    uint32_t md5::toLittleEndian32(uint32_t n) {
        if (!isBigEndian()) {
            return ((n << 24) & 0xFF000000) | ((n << 8) & 0x00FF0000) |
                ((n >> 8) & 0x0000FF00) | ((n >> 24) & 0x000000FF);
        }
        // Machine works on little endian, no need to change anything
        return n;
    }

    uint64_t md5::toLittleEndian64(uint64_t n) {
        if (!isBigEndian()) {
            return ((n << 56) & 0xFF00000000000000) |
                ((n << 40) & 0x00FF000000000000) |
                ((n << 24) & 0x0000FF0000000000) |
                ((n << 8) & 0x000000FF00000000) |
                ((n >> 8) & 0x00000000FF000000) |
                ((n >> 24) & 0x0000000000FF0000) |
                ((n >> 40) & 0x000000000000FF00) |
                ((n >> 56) & 0x00000000000000FF);
            ;
        }
        return n;
    }

    std::string md5::sig2hex(void* sig) {
        const char* hexChars = "0123456789abcdef";
        auto* intsig = static_cast<uint8_t*>(sig);
        std::string hex = "";
        for (uint8_t i = 0; i < 16; i++) {
            hex.push_back(hexChars[(intsig[i] >> 4) & 0xF]);
            hex.push_back(hexChars[(intsig[i]) & 0xF]);
        }
        return hex;
    }

    void* md5::hash_bs(const void* input_bs, uint64_t input_size) {
        auto* input = static_cast<const uint8_t*>(input_bs);

        // Step 0: Initial Data (Those are decided in the MD5 protocol)
        // s is the shift used in the leftrotate each round
        std::array<uint32_t, 64> s = {
            7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
            5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20,
            4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
            6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};
        // K is pseudo-random values used each round
        // The values can be obtained by the following python code:

        /**
        * @brief Values of K are pseudo-random and used to "salt" each round
        * The values can be obtained by the following python code
        * @code{.py}
        * from math import floor, sin
        *
        * for i in range(64):
        *     print(floor(2**32 * abs(sin(i+1))))
        * @endcode
        */
        std::array<uint32_t, 64> K = {
            3614090360, 3905402710, 606105819,  3250441966, 4118548399, 1200080426,
            2821735955, 4249261313, 1770035416, 2336552879, 4294925233, 2304563134,
            1804603682, 4254626195, 2792965006, 1236535329, 4129170786, 3225465664,
            643717713,  3921069994, 3593408605, 38016083,   3634488961, 3889429448,
            568446438,  3275163606, 4107603335, 1163531501, 2850285829, 4243563512,
            1735328473, 2368359562, 4294588738, 2272392833, 1839030562, 4259657740,
            2763975236, 1272893353, 4139469664, 3200236656, 681279174,  3936430074,
            3572445317, 76029189,   3654602809, 3873151461, 530742520,  3299628645,
            4096336452, 1126891415, 2878612391, 4237533241, 1700485571, 2399980690,
            4293915773, 2240044497, 1873313359, 4264355552, 2734768916, 1309151649,
            4149444226, 3174756917, 718787259,  3951481745};

        // The initial 128-bit state
        uint32_t a0 = 0x67452301, A = 0;
        uint32_t b0 = 0xefcdab89, B = 0;
        uint32_t c0 = 0x98badcfe, C = 0;
        uint32_t d0 = 0x10325476, D = 0;

        // Step 1: Processing the bytestring

        // First compute the size the padded message will have
        // so it is possible to allocate the right amount of memory
        uint64_t padded_message_size = 0;
        if (input_size % 64 < 56) {
            padded_message_size = input_size + 64 - (input_size % 64);
        } else {
            padded_message_size = input_size + 128 - (input_size % 64);
        }

        std::vector<uint8_t> padded_message(padded_message_size);

        // Beginning of the padded message is the original message
        std::copy(input, input + input_size, padded_message.begin());

        // Afterwards comes a single 1 bit and then only zeroes
        padded_message[input_size] = 1 << 7;  // 10000000
        for (uint64_t i = input_size; i % 64 != 56; i++) {
            if (i == input_size) {
                continue;  // pass first iteration
            }
            padded_message[i] = 0;
        }

        // We then have to add the 64-bit size of the message at the end
        // When there is a conversion from int to bytestring or vice-versa
        // We always need to make sure it is little endian
        uint64_t input_bitsize_le = toLittleEndian64(input_size * 8);
        for (uint8_t i = 0; i < 8; i++) {
            padded_message[padded_message_size - 8 + i] =
                (input_bitsize_le >> (56 - 8 * i)) & 0xFF;
        }

        // Already allocate memory for blocks
        std::array<uint32_t, 16> blocks{};

        // Rounds
        for (uint64_t chunk = 0; chunk * 64 < padded_message_size; chunk++) {
            // First, build the 16 32-bits blocks from the chunk
            for (uint8_t bid = 0; bid < 16; bid++) {
                blocks[bid] = 0;

                // Having to build a 32-bit word from 4-bit words
                // Add each and shift them to the left
                for (uint8_t cid = 0; cid < 4; cid++) {
                    blocks[bid] = (blocks[bid] << 8) +
                                padded_message[chunk * 64 + bid * 4 + cid];
                }
            }

            A = a0;
            B = b0;
            C = c0;
            D = d0;

            // Main "hashing" loop
            for (uint8_t i = 0; i < 64; i++) {
                uint32_t F = 0, g = 0;
                if (i < 16) {
                    F = (B & C) | ((~B) & D);
                    g = i;
                } else if (i < 32) {
                    F = (D & B) | ((~D) & C);
                    g = (5 * i + 1) % 16;
                } else if (i < 48) {
                    F = B ^ C ^ D;
                    g = (3 * i + 5) % 16;
                } else {
                    F = C ^ (B | (~D));
                    g = (7 * i) % 16;
                }

                // Update the accumulators
                F += A + K[i] + toLittleEndian32(blocks[g]);

                A = D;
                D = C;
                C = B;
                B += leftRotate32bits(F, s[i]);
            }
            // Update the state with this chunk's hash
            a0 += A;
            b0 += B;
            c0 += C;
            d0 += D;
        }

        // Build signature from state
        // Note, any type could be used for the signature
        // uint8_t was used to make the 16 bytes obvious
        // The sig needs to be little endian
        auto* sig = new uint8_t[16];
        for (uint8_t i = 0; i < 4; i++) {
            sig[i] = (a0 >> (8 * i)) & 0xFF;
            sig[i + 4] = (b0 >> (8 * i)) & 0xFF;
            sig[i + 8] = (c0 >> (8 * i)) & 0xFF;
            sig[i + 12] = (d0 >> (8 * i)) & 0xFF;
        }

        return sig;
    }

    void* md5::hash(const std::string& message) {
        return hash_bs(&message[0], message.size());
    }

}

namespace Tmdet::Utils {
    std::string Md5::getHash(std::string source) {
        auto md5 = hashing::md5();
        void* sig = md5.hash(source);
        return md5.sig2hex(sig);
    }
}