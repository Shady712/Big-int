#include "big_integer.h"
#include <cstddef>
#include <cstring>
#include <ostream>
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <functional>
#include <array>

constexpr int BASE = 10;
constexpr int BIT_SIZE = 32;

using int128 = __int128;

big_integer::big_integer() = default;

big_integer::big_integer(big_integer const& other) = default;

big_integer::big_integer(int a) : big_integer(static_cast<long long>(a)) {}


big_integer::big_integer(long a) : big_integer(static_cast<long long>(a)) {}

big_integer::big_integer(long long a) : minus(a < 0) {
    unsigned long long abs_a = a;
    while (abs_a > 0) {
        digits.push_back(abs_a);
        abs_a >>= BIT_SIZE;
    }
    delete_insignificant_digits();
}

big_integer::big_integer(unsigned int a) : big_integer(static_cast<unsigned long long>(a)) {}

big_integer::big_integer(unsigned long a) : big_integer(static_cast<unsigned long long>(a)) {}

big_integer::big_integer(unsigned long long a) {
    while (a > 0) {
        digits.push_back(a);
        a >>= BIT_SIZE;
    }
    delete_insignificant_digits();
}

bool is_digit(char x) {
    return '0' <= x && x <= '9';
}

uint32_t char_to_digit(char x) {
    return x - '0';
}

uint32_t get_number_from_substring(std::string const& str, size_t begin, size_t end) {
    assert(end <= str.size());
    uint32_t ans = 0;
    for (size_t i = begin; i < end; i++) {
        ans *= BASE;
        ans += char_to_digit(str[i]);
    }
    return ans;
}

big_integer::big_integer(std::string const& str) : big_integer() {
    if (str.empty() || (str.size() == 1 && !is_digit(str[0]))) {
        throw std::invalid_argument("Provided string is not a number");
    }
    for (size_t i = (str[0] == '+' || str[0] == '-') ? 1 : 0; i < str.size(); i++) {
        if (!is_digit(str[i])) {
            throw std::invalid_argument("Provided string is not a number");
        }
    }
    constexpr size_t SUBSTR_SIZE = 9;
    constexpr std::array<std::uint32_t , 10> powers{1, 10, 100, 1000, 10000,
                                                 100000, 1000000, 10000000, 100000000, 1000000000};
    for (size_t i = (str[0] == '+' || str[0] == '-') ? 1 : 0; i < str.size(); i += SUBSTR_SIZE) {
        *this *= powers[std::min(SUBSTR_SIZE, str.size() - i)];
        *this += get_number_from_substring(str, i, std::min(i + SUBSTR_SIZE, str.size()));
    }
    if (str[0] == '-') {
        negate();
    }
}

big_integer::~big_integer() = default;

uint32_t big_integer::get_insignificant_digit() const {
    return !minus ? 0 : UINT32_MAX;
}

uint32_t big_integer::get_digit(size_t index) const {
    return index < digits.size() ? digits[index] : get_insignificant_digit();
}

uint32_t get_negated_digit_with_carry(uint32_t digit, bool sign, uint32_t& carry) {
    if (!sign) {
        return digit;
    }
    digit ^= UINT32_MAX;
    digit += carry;
    if (digit > 0) {
        carry = 0;
    }
    return digit;
}

uint32_t from_64_to_32(uint64_t big_num) {
    big_num &= UINT32_MAX;
    return big_num;
}

void big_integer::add_index_with_carry(size_t index, uint32_t to_add, uint32_t& carry) {
    digits[index] += to_add;
    digits[index] += carry;
    if (carry > 0 && to_add == UINT32_MAX) {
        carry = 1;
    } else {
        if (digits[index] < to_add + carry) {
            carry = 1;
        } else {
            carry = 0;
        }
    }
}

bool big_integer::is_zero() const {
    return !minus && digits.empty();
}

bool big_integer::is_minus_one() const {
    return minus && digits.empty();
}

void big_integer::negate() {
    if (is_zero()) {
        return;
    }
    if (is_minus_one()) {
        minus = false;
        digits.push_back(1);
    } else {
        twos_complement_negate();
    }
}

void big_integer::invert_bits() {
    minus ^= true;
    for (uint32_t& i : digits) {
        i ^= UINT32_MAX;
    }
}

void big_integer::ensure_capacity(size_t cap) {
    if (digits.size() < cap) {
        digits.resize(cap, get_insignificant_digit());
    }
}

void big_integer::delete_insignificant_digits() {
    uint32_t true_zero = get_insignificant_digit();
    while (!digits.empty() && digits.back() == true_zero) {
        digits.pop_back();
    }
}

void big_integer::twos_complement_negate() {
    invert_bits();
    *this += 1;
    delete_insignificant_digits();
}

big_integer& big_integer::operator=(big_integer const& other) = default;

big_integer& big_integer::operator+=(big_integer const& rhs) {
    if (rhs.is_zero()) {
        return *this;
    }

    ensure_capacity(std::max(digits.size(), rhs.digits.size()) + 1);
    uint32_t carry = 0;
    for (size_t i = 0; i < digits.size(); i++) {
        uint32_t digit_to_add = rhs.get_digit(i);
        add_index_with_carry(i, digit_to_add, carry);
    }
    minus = (digits.back() & (1 << (BIT_SIZE - 1))) > 0;
    delete_insignificant_digits();
    return *this;
}

big_integer& big_integer::operator-=(big_integer const& rhs) {
    if (rhs.is_zero()) {
        return *this;
    }
    if (rhs.is_minus_one()) {
        return *this += 1;
    }

    ensure_capacity(std::max(digits.size(), rhs.digits.size()) + 1);
    uint32_t carry = 0;
    uint32_t negate_carry = 1;
    for (size_t i = 0; i < digits.size(); i++) {
        uint32_t digit_to_add = get_negated_digit_with_carry(rhs.get_digit(i), true, negate_carry);
        add_index_with_carry(i, digit_to_add, carry);
    }

    if (negate_carry > 0) {
        uint32_t digit_to_add = get_negated_digit_with_carry(rhs.get_insignificant_digit(), true, negate_carry);
        add_index_with_carry(digits.size() - 1, digit_to_add, carry);
    }

    minus = (digits.back() & (1 << (BIT_SIZE - 1))) > 0;
    delete_insignificant_digits();
    return *this;
}

void big_integer::mul_iteration(uint64_t digit_a, const big_integer& b, big_integer& ans, size_t i) {
    uint32_t digit_ans = 0;
    uint32_t negate_carry_b = b.minus;
    for (size_t j = 0; j < b.digits.size() + negate_carry_b; j++) {
        uint64_t digit_b = get_negated_digit_with_carry(b.get_digit(j), b.minus, negate_carry_b);
        uint64_t addition = ans.get_digit(i + j);
        uint64_t product = digit_a * digit_b + addition + digit_ans;
        digit_ans = product >> BIT_SIZE;
        ans.digits[i + j] = from_64_to_32(product);
    }
    ans.digits[i + b.digits.size()] = digit_ans;
}

big_integer& big_integer::operator*=(big_integer const& rhs) {
    big_integer ans;
    ans.ensure_capacity(digits.size() + rhs.digits.size() + 1);
    uint32_t negate_carry = minus;
    for (size_t i = 0; i < digits.size() + negate_carry; i++) {
        mul_iteration(get_negated_digit_with_carry(get_digit(i), minus, negate_carry), rhs, ans, i);
    }
    ans.delete_insignificant_digits();
    if (minus != rhs.minus) {
        ans.negate();
    }
    std::swap(*this, ans);
    return *this;
}

uint32_t big_integer::short_divide(uint32_t x) {
    bool negated = minus;
    if (negated) {
        negate();
    }
    ensure_capacity(digits.size());
    const size_t end = digits.size() - 1;
    uint64_t carry = 0;
    for (size_t i = 0; i < digits.size(); i++) {
        uint64_t cur = digits[end - i] | (carry << BIT_SIZE);
        digits[end - i] = from_64_to_32(cur / x);
        carry = cur % x;
    }
    delete_insignificant_digits();
    if (negated) {
        negate();
    }
    return carry;
}

bool big_integer::more_than_shifted(big_integer const& other, size_t shift) {
    return !less_than_shifted(other, shift);
}

bool big_integer::less_than_shifted(big_integer const& other, size_t shift) {
    for (size_t i = digits.size(); i > 0; i--) {
        size_t shift_id = shift - digits.size() + i - 1;
        if (other.get_digit(shift_id) != digits[i - 1]) {
            return other.get_digit(shift_id) < digits[i - 1];
        }
    }
    return true;
}

void big_integer::shifted_subtract(big_integer const& other, size_t shift) {
    uint32_t carry = 0;
    uint32_t negate_carry = 1;
    for (size_t i = 0; i < shift; i++) {
        uint32_t digit_to_add = get_negated_digit_with_carry(other.get_digit(i), true, negate_carry);
        add_index_with_carry(digits.size() - shift + i, digit_to_add, carry);
    }
    if (digits.back() == 0) {
        digits.pop_back();
    }
}

int128 big_integer::get_shifted_digit(size_t bit_count) {
    int128 ans = 0;
    for (size_t j = 1; j <= bit_count; j++) {
        int128 cur = get_digit(digits.size() - j);
        cur <<= BIT_SIZE * (bit_count - j);
        ans |= cur;
    }
    return ans;
}

big_integer& big_integer::operator/=(big_integer const& rhs) {
    bool true_minus = minus;
    if (minus) {
        negate();
    }
    big_integer abs_b = rhs.minus ? -rhs : rhs;
    if (*this < abs_b) {
        return *this = 0;
    }
    if (abs_b.digits.size() == 1) {
        short_divide(abs_b.get_digit(0));
        if (true_minus != rhs.minus) {
            negate();
        }
        return *this;
    }
    const size_t size_a = digits.size() + 1;
    const size_t size_b = abs_b.digits.size() + 1;
    digits.push_back(0);
    big_integer ans;
    ans.ensure_capacity(size_a - size_b + 1);
    size_t index_ans = size_a - size_b;
    for (size_t i = size_b; i <= size_a; i++) {
        int128 digit_a = get_shifted_digit(3);
        int128 digit_b = abs_b.get_shifted_digit(2);
        uint32_t digit_tmp = static_cast<uint32_t>(std::min(static_cast<int128>(UINT32_MAX), digit_a / digit_b));
        big_integer tmp = abs_b * digit_tmp;
        if (more_than_shifted(tmp, size_b)) {
            tmp -= abs_b;
            digit_tmp--;
        }
        shifted_subtract(tmp, size_b);
        ans.digits[index_ans--] = digit_tmp;
    }
    ans.delete_insignificant_digits();
    if (true_minus != rhs.minus) {
        ans.negate();
    }
    std::swap(*this, ans);
    return *this;
}

big_integer& big_integer::operator%=(big_integer const& rhs) {
    *this = *this - *this / rhs * rhs;
    return *this;
}

uint32_t bit_and(uint32_t x, uint32_t y) {
    return x & y;
}

uint32_t bit_or(uint32_t x, uint32_t y) {
    return x | y;
}

uint32_t bit_xor(uint32_t x, uint32_t y) {
    return x ^ y;
}

big_integer& big_integer::abstract_bit_operation(big_integer const& rhs, const std::function<uint32_t(uint32_t, uint32_t)>& operation) {
    ensure_capacity(rhs.digits.size());
    for (size_t i = 0; i < digits.size(); i++) {
        digits[i] = operation(digits[i], rhs.get_digit(i));
    }
    minus = operation(minus, rhs.minus);
    delete_insignificant_digits();
    return *this;
}

big_integer& big_integer::operator&=(big_integer const& rhs) {
    return abstract_bit_operation(rhs, bit_and);
}

big_integer& big_integer::operator|=(big_integer const& rhs) {
    return abstract_bit_operation(rhs, bit_or);
}

big_integer& big_integer::operator^=(big_integer const& rhs) {
    return abstract_bit_operation(rhs, bit_xor);
}

big_integer& big_integer::operator<<=(int rhs) {
    size_t new_digits = rhs / BIT_SIZE;
    ensure_capacity(digits.size() + new_digits);
    for (size_t i = digits.size(); i > new_digits; i--) {
        size_t index = i - 1;
        digits[index] = get_digit(index - new_digits);
    }
    for (size_t i = 0; i < new_digits; i++) {
        digits[i] = 0;
    }
    delete_insignificant_digits();
    if (rhs % BIT_SIZE > 0) {
        digits.insert(digits.begin(), 0);
        *this >>= BIT_SIZE - rhs % BIT_SIZE;
    }
    return *this;
}

big_integer& big_integer::operator>>=(int rhs) {
    size_t deleted_digits = rhs / BIT_SIZE;
    for (size_t i = deleted_digits; i < digits.size(); i++) {
        digits[i - deleted_digits] = get_digit(i);
    }
    for (size_t i = 0; i < deleted_digits; i++) {
        digits.pop_back();
    }
    delete_insignificant_digits();
    uint32_t remainder = short_divide((static_cast<uint32_t>(1)) << (rhs % BIT_SIZE));
    if (*this < 0 && remainder > 0) {
        *this -= 1;
    }
    return *this;
}

big_integer big_integer::operator+() const {
    return *this;
}

big_integer big_integer::operator-() const {
    big_integer ans = *this;
    ans.negate();
    return ans;
}

big_integer big_integer::operator~() const {
    big_integer ans = *this;
    ans.invert_bits();
    return ans;
}

big_integer& big_integer::operator++() {
    return *this += 1;
}

big_integer big_integer::operator++(int) {
    big_integer ans = *this;
    *this += 1;
    return ans;
}

big_integer& big_integer::operator--() {
    return *this -= 1;
}

big_integer big_integer::operator--(int) {
    big_integer ans = *this;
    *this -= 1;
    return ans;
}

big_integer operator+(big_integer a, big_integer const& b) {
    return a += b;
}

big_integer operator-(big_integer a, big_integer const& b) {
    return a -= b;
}

big_integer operator*(big_integer a, big_integer const& b) {
    return a *= b;
}

big_integer operator/(big_integer a, big_integer const& b) {
    return a /= b;
}

big_integer operator%(big_integer a, big_integer const& b) {
    return a %= b;
}

big_integer operator&(big_integer a, big_integer const& b) {
    return a &= b;
}

big_integer operator|(big_integer a, big_integer const& b) {
    return a |= b;
}

big_integer operator^(big_integer a, big_integer const& b) {
    return a ^= b;
}

big_integer operator<<(big_integer a, int b) {
    return a <<= b;
}

big_integer operator>>(big_integer a, int b) {
    return a >>= b;
}

bool operator==(big_integer const& a, big_integer const& b) {
    return a.minus == b.minus && a.digits == b.digits;
}

bool operator!=(big_integer const& a, big_integer const& b) {
    return !(a == b);
}

bool operator<(big_integer const& a, big_integer const& b) {
    if (a.minus == b.minus) {
        if (a.digits.size() == b.digits.size()) {
            for (size_t i = a.digits.size(); i > 0; i--) {
                size_t id = i - 1;
                if (a.digits[id] != b.digits[id]) {
                    return a.digits[id] < b.digits[id];
                }
            }
            return false;
        } else {
            if (a.minus) {
                return a.digits.size() > b.digits.size();
            }
            else {
                return a.digits.size() < b.digits.size();
            }
        }
    } else {
        return a.minus;
    }
}

bool operator>(big_integer const& a, big_integer const& b) {
    return b < a;
}

bool operator<=(big_integer const& a, big_integer const& b) {
    return !(b < a);
}

bool operator>=(big_integer const& a, big_integer const& b) {
    return !(a < b);
}

char digit_to_char(uint32_t digit) {
    return '0' + digit;
}

std::string to_string(big_integer const& a) {
    if (a == 0) {
        return "0";
    }
    big_integer abs_a = a.minus ? -a : a;
    std::string ans;
    constexpr int DIGIT_BASE = 1000000000;
    while (abs_a > 0) {
        uint32_t digit = abs_a.short_divide(DIGIT_BASE);
        for (size_t i = 0; i < BASE - 1; i++) {
            ans += digit_to_char(digit % BASE);
            digit /= BASE;
        }
    }
    while (!ans.empty() && ans.back() == '0') {
        ans.pop_back();
    }
    if (a.minus) {
        ans += '-';
    }
    std::reverse(ans.begin(), ans.end());
    return ans;
}

std::ostream& operator<<(std::ostream& s, big_integer const& a) {
    return s << to_string(a);
}
