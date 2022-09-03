#include <iostream>
#include <queue>
#include <bits/stdc++.h>
#include <vector>

using namespace std;

vector<int> x;

vector<int> testando()
{
    x.push_back(1);    
}

int main(int argc, char *argv[])
{
    auto b = testando();
    cout << "b " << b.size() << endl;
    b.clear();
    cout << "x " << x.size() << endl;
    x.clear();

    return 0;
}
