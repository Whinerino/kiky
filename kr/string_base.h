#ifndef STRING_BASE_H
#define STRING_BASE_H

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>

// Базовый класс
class String {
protected:
    std::string data;
public:
    // Конструкторы
    String();
    String(const std::string& str);
    
    // Виртуальный деструктор
    virtual ~String();

    // Чисто виртуальная функция перегрузки оператора +=
    virtual String& operator+=(char c) = 0;

    // Перегрузка оператора вывода
    friend std::ostream& operator<<(std::ostream& os, const String& s);
};

// Дочерний класс String1 (по 1 символу в начало и конец)
class String1 : public String {
public:
    String1(const std::string& str);
    String& operator+=(char c);  // ← Убрали override
};

// Дочерний класс String2 (по 2 символа в начало и конец)
class String2 : public String {
public:
    String2(const std::string& str);
    String& operator+=(char c);  // ← Убрали override
};

// Вспомогательная функция для проверки, состоит ли слово только из цифр
bool isAllDigits(const std::string& str);

#endif // STRING_BASE_H