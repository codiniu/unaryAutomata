This code allows to generate all non-isomorphic unary automata.

To compile the program with default options, just call: **make**.  
To generate unary automata of N states, call generate() function on unaryAutomata<N> class object.  
Number of states has to be known at compile time.  
To print the output, call **./unaryAutomata -p**.

In constants.hpp you can also set available flags:
- maxNrOfCycles - generate automata with number of cycles less or equal to maxNrOfCycles
- minNrOfCycles - generate automata with number of cycles greater or equal to minNrOfCycles

Dawid Wójcik (wojcikdawid92@gmail.com)
