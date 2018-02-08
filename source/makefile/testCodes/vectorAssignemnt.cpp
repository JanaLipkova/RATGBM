// vector assignment
#include <iostream>
#include <vector>
#include <stdio.h>
#include <assert.h>     /* assert */

int main ()
{
  std::vector<int> foo (3,0);
  std::vector<int>* pfoo ;
  std::vector<int> bar ;
  std::vector<int> tmp = ( 1, 2, 3);
  //bar = foo;
  //foo = std::vector<int>();
  


  std::cout << "Size of foo: " << int(foo.size()) << '\n';
  std::cout << "Size of bar: " << int(bar.size()) << '\n';
  

 foo[0] = 1;
 foo[1] = 2;
 foo[2] = 3;

  pfoo = &foo;
  bar = *pfoo;

  //bar = foo;

  std::cout << "Size of foo: " << int(foo.size()) << '\n';
  std::cout << "Size of bar: " << int(bar.size()) << '\n';

 printf("foo=[%d,%d,%d], \n ", foo[0], foo[1], foo[2]);
 printf("bar=[%d,%d,%d], \n ", bar[0], bar[1], bar[2]);
  
 assert(1 && 1); 
 
 return 0;
}
