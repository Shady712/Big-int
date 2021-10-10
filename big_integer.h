#pragma once

#include <iosfwd>
#include <string>
#include <vector>
#include <functional>

using int128 = __int128;

struct big_integer {
    big_integer();
    big_integer(big_integer const& other);
    big_integer(unsigned long long a);
    big_integer(unsigned long a);
    big_integer(unsigned int a);
    big_integer(long long a);
    big_integer(long a);
    big_integer(int a);
    explicit big_integer(std::string const& str);
    ~big_integer();

    big_integer& operator=(big_integer const& other);

    big_integer& operator+=(big_integer const& rhs);
    big_integer& operator-=(big_integer const& rhs);
    big_integer& operator*=(big_integer const& rhs);
    big_integer& operator/=(big_integer const& rhs);
    big_integer& operator%=(big_integer const& rhs);

    big_integer& operator&=(big_integer const& rhs);
    big_integer& operator|=(big_integer const& rhs);
    big_integer& operator^=(big_integer const& rhs);

    big_integer& operator<<=(int rhs);
    big_integer& operator>>=(int rhs);

    big_integer operator+() const;
    big_integer operator-() const;
    big_integer operator~() const;

    big_integer& operator++();
    big_integer operator++(int);

    big_integer& operator--();
    big_integer operator--(int);


    friend big_integer operator+(big_integer a, big_integer const& b);
    friend big_integer operator-(big_integer a, big_integer const& b);
    friend big_integer operator*(big_integer a, big_integer const& b);
    friend big_integer operator/(big_integer a, big_integer const& b);
    friend big_integer operator%(big_integer a, big_integer const& b);

    friend big_integer operator<<(big_integer a, int b);
    friend big_integer operator>>(big_integer a, int b);

    friend bool operator==(big_integer const& a, big_integer const& b);
    friend bool operator!=(big_integer const& a, big_integer const& b);
    friend bool operator<(big_integer const& a, big_integer const& b);
    friend bool operator>(big_integer const& a, big_integer const& b);
    friend bool operator<=(big_integer const& a, big_integer const& b);
    friend bool operator>=(big_integer const& a, big_integer const& b);

    friend std::string to_string(big_integer const& a);

private:
    bool minus = false;
    std::vector<uint32_t> digits;

    void negate();
    void invert_bits();
    void twos_complement_negate();

    uint32_t get_insignificant_digit() const;
    uint32_t get_digit(size_t index) const;
    void delete_insignificant_digits();

    void add_index_with_carry(size_t index, uint32_t to_add, uint32_t& carry);

    bool is_zero() const;
    bool is_minus_one() const;

    void ensure_capacity(size_t cap);

    static void mul_iteration(uint64_t digit_a, big_integer const& b, big_integer& ans, size_t i);

    uint32_t short_divide(uint32_t x);
    bool more_than_shifted(big_integer const& other, size_t shift);
    bool less_than_shifted(big_integer const& other, size_t shift);
    void shifted_subtract(big_integer const& other, size_t shift);
    int128 get_shifted_digit(size_t bit_count);

    big_integer& abstract_bit_operation(big_integer const& rhs, std::function<uint32_t(uint32_t, uint32_t)> const& operation);
};


big_integer operator&(big_integer a, big_integer const& b);
big_integer operator|(big_integer a, big_integer const& b);
big_integer operator^(big_integer a, big_integer const& b);

bool operator==(big_integer const& a, big_integer const& b);
bool operator!=(big_integer const& a, big_integer const& b);
bool operator<(big_integer const& a, big_integer const& b);
bool operator>(big_integer const& a, big_integer const& b);
bool operator<=(big_integer const& a, big_integer const& b);
bool operator>=(big_integer const& a, big_integer const& b);

std::string to_string(big_integer const& a);
std::ostream& operator<<(std::ostream& s, big_integer const& a);
