#include <iostream>
#include <cstdlib>
#include <string>

#define ASSERT_EQUAL(val, expected) \
    do { \
        auto _val = (val); \
        auto _expected = (expected); \
        std::string _val_str = std::to_string(_val); \
        std::string _expected_str = std::to_string(_expected); \
        if (_val_str == _expected_str) { \
            std::cout << "PASS: " << #val << " == " << #expected << " (" << _val_str << ")\n"; \
        } else { \
            std::cerr << "FAIL: " << #val << " != " << #expected << " (" << _val_str << " != " << _expected_str << ")\n"; \
            size_t len = std::min(_val_str.size(), _expected_str.size()); \
            size_t match = 0; \
            while (match < len && _val_str[match] == _expected_str[match]) ++match; \
            std::cerr << "First difference at digit position: " << match << "\n"; \
            std::cerr << "  " << #val << ":      " << _val_str << "\n"; \
            std::cerr << "  " << #expected << ": " << _expected_str << "\n"; \
            std::cerr << "  "; \
            for (size_t i = 0; i < match; ++i) std::cerr << " "; \
            std::cerr << "^\n"; \
            std::abort(); \
        } \
    } while (0)
