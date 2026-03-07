#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "string_base.h"

int main() {
    // Открываем файл для чтения
    std::ifstream file("data.dat");
    if (!file.is_open()) {
        std::cerr << "Ошибка: Не удалось открыть файл data.dat" << std::endl;
        return 1;
    }

    std::vector<String*> strings;
    std::string word;

    // Читаем файл слово за словом (автоматически игнорируются пробелы, \t, \n)
    while (file >> word) {
        String* newString = nullptr;

        // Если только цифры — создаем String1, иначе String2
        if (isAllDigits(word)) {
            newString = new String1(word);
        } else {
            newString = new String2(word);
        }

        // Применяем оператор += с символом '\' (обратный слэш)
        *newString += '\\';
        
        // Добавляем указатель в вектор
        strings.push_back(newString);
    }
    file.close();

    // Выводим результаты на экран и освобождаем память
    std::cout << "Результаты:" << std::endl;
    for (String* s : strings) {
        std::cout << *s << std::endl;
        delete s;  // Важно: предотвратить утечку памяти
    }

    return 0;
}