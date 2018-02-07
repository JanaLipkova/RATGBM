#include <iostream>
int main()
{
	enum Color { RED, GREEN, BLUE };
	Color r = RED;

switch(r)
{
    case RED  : std::cout << "red\n";   break;
    case GREEN: std::cout << "green\n"; break;
    case BLUE : std::cout << "blue\n";  break;
}

}
