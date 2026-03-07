#include "string_base.h"

// ===== Вспомогательная функция (вынесена отдельно для C++98) =====
// Функция для проверки, является ли символ цифрой
static bool isDigitChar(unsigned char c) {
    return std::isdigit(c) != 0;
}

// ===== Реализация базового класса String =====

String::String() : data("") {}

String::String(const std::string& str) : data(str) {}

String::~String() {}

// Оператор вывода (определение вне класса)
std::ostream& operator<<(std::ostream& os, const String& s) {
    os << s.data;
    return os;
}

// ===== Реализация класса String1 =====

String1::String1(const std::string& str) : String(str) {}

String& String1::operator+=(char c) {
    data = c + data + c;
    return *this;
}

// ===== Реализация класса String2 =====

String2::String2(const std::string& str) : String(str) {}

String& String2::operator+=(char c) {
    // Создаем строку из 2-х символов c (совместимо с C++98)
    std::string append_str;
    append_str += c;
    append_str += c;
    data = append_str + data + append_str;
    return *this;
}

// ===== Вспомогательная функция =====

bool isAllDigits(const std::string& str) {
    if (str.empty()) return false;
    
    // Используем std::find_if вместо лямбды (C++98 совместимо)
    return std::find_if(str.begin(), str.end(), isDigitChar) == str.end();
}