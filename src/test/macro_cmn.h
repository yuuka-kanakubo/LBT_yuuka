#define ASSERT_EQUAL(val, expected) \
    do { \
        if ((val) == (expected)) { \
            std::cout << "PASS: " << #val << " == " << #expected << " (" << val << ")\n"; \
        } else { \
            std::cerr << "FAIL: " << #val << " != " << #expected << " (" << val << " != " << expected << ")\n"; \
            std::abort(); \
        } \
    } while (0)
