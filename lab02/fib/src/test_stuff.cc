#include <gtest/gtest.h>
#include "matrix.h"
#include "compute.cpp"
#include "compute2.cpp"
#include "computeNew.cpp"

TEST(MatrixVectorTest, ComputeTest) {
    std::vector<int> v = {10, 2, 20};
    Matrix m(3, 3);
    
    for (int i = 0; i < 3; i++) {
        for (int y = 0; y < 3; y++) {
            m(i, y) = 1;
        }
    }

    Matrix expected(3, 3);
    for (int i = 0; i < 3; i++) {
        expected(0, i) = 31.7662;
    }
    for (int i = 0; i < 3; i++) {
        expected(1, i) = 0.82224;
    }
    for (int i = 0; i < 3; i++) {
        expected(2, i) = -1.6214;
    }
    
    auto erg = compute(m, v);
    
    // Erwartete Ergebnisse für die Multiplikation von Matrix und Vektor
    EXPECT_EQ(v[0], 10);
    EXPECT_EQ(v[1], 2);
    EXPECT_EQ(v[2], 20);
    
    // Testen, ob das Ergebnis der Matrix korrekt ist
    for (int i = 0; i < 3; i++) {
        for (int y = 0; y < 3; y++) {
            EXPECT_NEAR(erg(i,y), expected(i, y), 0.1); // Falls es sich um eine einfache Summation handelt
        }
    }
}

TEST(MatrixVectorTest, Compute2Test) {
    std::vector<int> v = {10, 2, 20};
    Matrix m(3, 3);
    
    for (int i = 0; i < 3; i++) {
        for (int y = 0; y < 3; y++) {
            m(i, y) = 1;
        }
    }

    Matrix expected(3, 3);
    for (int i = 0; i < 3; i++) {
        expected(0, i) = 31.7662;
    }
    for (int i = 0; i < 3; i++) {
        expected(1, i) = 0.82224;
    }
    for (int i = 0; i < 3; i++) {
        expected(2, i) = -1.6214;
    }
    
    auto erg = compute2(m, v);
    
    // Erwartete Ergebnisse für die Multiplikation von Matrix und Vektor
    EXPECT_EQ(v[0], 10);
    EXPECT_EQ(v[1], 2);
    EXPECT_EQ(v[2], 20);
    
    // Testen, ob das Ergebnis der Matrix korrekt ist
    for (int i = 0; i < 3; i++) {
        for (int y = 0; y < 3; y++) {
            EXPECT_NEAR(erg(i,y), expected(i, y), 0.1); // Falls es sich um eine einfache Summation handelt
        }
    }
}

TEST(MatrixVectorTest, ComputeNewTest) {
    std::vector<int> v = {10, 2, 20};
    Matrix m(3, 3);
    
    for (int i = 0; i < 3; i++) {
        for (int y = 0; y < 3; y++) {
            m(i, y) = 1;
        }
    }

    Matrix expected(3, 3);
    for (int i = 0; i < 3; i++) {
        expected(0, i) = 31.7662;
    }
    for (int i = 0; i < 3; i++) {
        expected(1, i) = 0.82224;
    }
    for (int i = 0; i < 3; i++) {
        expected(2, i) = -1.6214;
    }
    
    auto erg = computeNew(m, v);
    
    // Erwartete Ergebnisse für die Multiplikation von Matrix und Vektor
    EXPECT_EQ(v[0], 10);
    EXPECT_EQ(v[1], 2);
    EXPECT_EQ(v[2], 20);
    
    // Testen, ob das Ergebnis der Matrix korrekt ist
    for (int i = 0; i < 3; i++) {
        for (int y = 0; y < 3; y++) {
            EXPECT_NEAR(erg(i,y), expected(i, y), 0.1); // Falls es sich um eine einfache Summation handelt
        }
    }
}
