#include "unaryAutomata.hpp"

int main(int argc, char *argv[])
{
  if (argc == 2)
  {
    std::string arg{argv[1]};
    if (arg == "-p")
      printAutomata = true;
  }

  unaryAutomata<8> ua;
  ua.generate();

  return 0;
}