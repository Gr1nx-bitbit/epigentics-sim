#include <iostream>
#include <string>
using namespace std;

int main(void) {
    string hello = "hey";
    if (hello[3] == '\0') {
        cout << "its true!" << endl;
    }
    return 0;
}