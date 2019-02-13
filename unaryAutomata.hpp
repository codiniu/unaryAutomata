#ifndef UNARYAUTOMATA_HPP_
#define UNARYAUTOMATA_HPP_

#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <numeric>

#include "constants.hpp"
#include "trees.hpp"

void reportAutomaton(uint8_t* aut, uint size)
{
  if (printAutomata)
  {
    for (uint i = 0; i < size; i++)
      std::cout << +aut[i] << " ";
    std::cout << std::endl;
  }
}

template <uint N>
class unaryAutomata
{
  TreesUpToN<N/2> trees;
  std::vector<uint> combs[N/2+1];
  uint64_t autCnt;
  uint64_t automataNr = 0;
  uint automataSize = 0;
  uint tailCnt = 0;
  uint8_t* automata;
  uint8_t* automataTail;
  bool hasMoreThanHalfVs = false;

public:
  unaryAutomata()
  : autCnt(treesNr[N]+1),
  automata(new uint8_t[connAutUpToN[N/2-1] * (N/2)]),
  automataTail(new uint8_t[unaryAutomataNr[N/2]*N/2])
  {
    if (minNrOfCycles <= 1)
    {
      TreesIterator treesIt(N);
      uint8_t* currTree;
      while ((currTree = treesIt.next()))
        reportAutomaton(currTree, N);
      autCnt = treesNr[N]+1;
    }

    else
      autCnt = 1;

    for (uint i=1; i <= N/2; i++)
    {
      automataSize++;
      for (uint j=1; j < i; j++)
      {
        generatePartitions(i-j, j);
      }

      uint8_t lastAutomaton[N]{};
      for (uint k=0; k < automataSize; k++)
        lastAutomaton[k] = (k+1)%i;
      std::memcpy((*this)[automataNr], lastAutomaton, N/2*sizeof(uint8_t));
      automataNr++;
    }
  }

  inline uint8_t* operator[](const uint64_t k) { return automata + N/2 * k; }
  inline const uint8_t* operator[](const uint64_t k) const { return automata + N/2 * k; }

  void reportAutomatonWithTails(uint8_t* aut)
  {
    uint tailSize = N - automataSize;
    autCnt += unaryAutomataNr[tailSize];
    if (printAutomata)
    {
      for (uint i = 0; i < unaryAutomataNr[tailSize]; i++)
      {
        for (uint j = 0; j < automataSize; j++)
        {
          std::cout << +aut[j] << " ";
        }
        for (uint k = 0; k < tailSize; k++)
          std::cout << +automataTail[i*(N/2) + k] << " ";
        std::cout << std::endl;
      }
    }
  }

  uint getPeriodLength(uint* a, uint n) 
  {
    uint lps[N];
    uint len = 0;
    uint i = 1;

    lps[0] = 0;
    while (i < n) 
    {
      if (a[i] == a[len])
        lps[i++] = ++len; 
      else
      {
        if (len != 0) 
          len = lps[len-1];
        else
          lps[i++] = 0; 
      }
    }
    len = lps[n-1]; 
    return n-len;
  }

  void generatePartitions(const uint n, const uint l)
  {
    if (n == N-1)
      return;

    uint a[N]{};
    uint k = 1;
    a[1] = n;
    while (k != 0)
    {
      uint x = a[k-1] + 1;
      uint y = a[k] - 1;
      k--;
      while (x <= y && k < l-1)
      {
        a[k] = x;
        y -= x;
        k++;
      }
      a[k] = x + y;
      generateNecklaces(a, k, l);
    }
  }

  void generateNecklaces(const uint* a, const uint end, const uint cycleSize)
  {
    uint result[N]{};
    uint M[N]{};
    uint c[N]{};
    uint cnt[N]{};
    uint k = 0;
    uint n = 0;

    cnt[0] = cycleSize-end-1;
    for (uint i=0; i <= end; i++)
      cnt[a[i]]++;
    for (uint i=0; i < automataSize; i++)
    {
      if (cnt[i] != 0)
      {
        M[k] = i;
        c[k] = cnt[i];
        n += cnt[i];
        k++;
      }
    }
    c[0]--;
    necklacesWithFixedContent(2, 1, result, M, c, k, n);
  }

  void necklacesWithFixedContent(uint t, uint p, uint* res, uint* M, uint* c,
    uint k, uint n)
  {
    if (t > n && n%p == 0)
    {
      if (n!=p && n > 1 && M[k-1] > 1)
        generateUnaryAutomata(M, res, n, true);
      else
        generateUnaryAutomata(M, res, n, false);
    }
    else
    {
      for (uint j = res[t-p]; j < k; j++)
      {
        if (c[j])
        {
          res[t] = j;
          c[j]--;
          if (j == res[t-p])
            necklacesWithFixedContent(t+1, p, res, M, c, k, n);
          else
            necklacesWithFixedContent(t+1, t, res, M, c, k, n);
          c[j]++;
        }
      }
    }
  }

  void generateUnaryAutomata(uint* M, uint* necklace, uint n, bool isCyclic)
  {
    uint8_t automaton[N]{};
    uint cycle[N]{};
    uint cnt[N]{};
    uint cyclePos = 0;

    for (uint i = 1; i <= n; i++)
    {
      uint nrOfVertices = M[necklace[i]];
      cnt[nrOfVertices]++;
      uint nextPos = cyclePos + nrOfVertices + 1;
      if (nextPos == automataSize)
        automaton[cyclePos] = 0;
      else
      {
        automaton[cyclePos] = nextPos;
        cycle[i] = nextPos;
        cyclePos += nrOfVertices + 1;
      }
    }

    uint mappedNecklace[N]{};
    for (uint i = 0; i < n; i++)
      mappedNecklace[i] = M[necklace[i+1]];

    if (isCyclic)
    {
      TrieNode* root = new TrieNode(treesNr[mappedNecklace[0]+1]);
      uint currConf[N]{};
      uint periodLength = getPeriodLength(mappedNecklace, n);
      generateAutomataWithCyclicNecklace(0, automaton, mappedNecklace, cycle, periodLength, n, currConf, root);
      delete root;
    }
    else
    {
      attachTreesToCycle(0, automaton, mappedNecklace, cycle, n);
    }
  }

  void generateAutomataWithCyclicNecklace(uint level, uint8_t* automaton, uint* necklace, uint* cycle, uint periodLength, uint n,
    uint* currConf, TrieNode* chosen)
  {
    const uint nrOfVertices = necklace[level];

    if (n == 2)
    {
      uint start = treesNrUpToN[nrOfVertices];
      uint end = treesNrUpToN[nrOfVertices+1];

      for (uint i = start; i < end; i++)
      {
        currConf[0] = i;
        for (uint j = i; j < end; j++)
        {
          currConf[1] = j;
          for (uint k = 0; k < n; k++)
          {
            for (uint l = 1; l <= necklace[0]; l++)
            {
              uint cyclePos = cycle[k];
              automaton[cyclePos+l] = trees[currConf[k]][l] + cyclePos;
            }
          }
          if (!hasMoreThanHalfVs)
          {
            std::memcpy((*this)[automataNr], automaton, N/2*sizeof(uint8_t));
            automataNr++;
          }
          else
          {
            if (automataSize == N){
              reportAutomaton(automaton, automataSize);
              autCnt++;
            }
            else
              reportAutomatonWithTails(automaton);
          }
        }
      }
      return;
    }

    if (level == n)
    {
      if (!isConfInTrie(currConf, necklace, n, chosen))
      {
        insertConf(currConf, necklace, n, chosen);
        for (uint j = 0; j < n/periodLength - 1; j++)
        {
          uint rotation[N];
          for (uint k = 0; k < n; k++)
            rotation[k] = currConf[(k + n - (j+1)*periodLength)%n];
          insertConf(rotation, necklace, n, chosen);
        }

        for (uint i = 0; i < n; i++)
        {
          for (uint j = 1; j <= necklace[i]; j++)
          {
            uint cyclePos = cycle[i];
            automaton[cyclePos+j] = trees[currConf[i]][j] + cyclePos;
          }
        }

        if (!hasMoreThanHalfVs)
        {
          std::memcpy((*this)[automataNr], automaton, N/2*sizeof(uint8_t));
          automataNr++;
        }
        else
        {
          if (automataSize == N){
            reportAutomaton(automaton, automataSize);
            autCnt++;
          }
          else
            reportAutomatonWithTails(automaton);
        }
      }
      return;
    }

    for (uint i = treesNrUpToN[nrOfVertices]; i < treesNrUpToN[nrOfVertices+1]; i++)
    {
      currConf[level] = i;
      generateAutomataWithCyclicNecklace(level+1, automaton, necklace, cycle, periodLength, n, currConf, chosen);
    } 
  }

  void attachTreesToCycle(uint level, uint8_t* automaton, uint* necklace, uint* cycle, uint n)
  {
    if (level == n)
    {
      if (!hasMoreThanHalfVs)
      {
        std::memcpy((*this)[automataNr], automaton, N/2*sizeof(uint8_t));
        automataNr++;
      }
      else
      {
        if (automataSize == N){
          reportAutomaton(automaton, automataSize);
          autCnt++;
        }
        else
          reportAutomatonWithTails(automaton);
      }
      return;
    }

    uint cyclePos = cycle[level];
    const uint nrOfVertices = necklace[level];

    if (nrOfVertices < N/2)
    {
      for (uint i = treesNrUpToN[nrOfVertices]; i < treesNrUpToN[nrOfVertices+1]; i++)
      {
        for (uint j = 1; j <= nrOfVertices; j++)
        {
          automaton[cyclePos+j] = trees[i][j] + cyclePos;
        }
        attachTreesToCycle(level+1, automaton, necklace, cycle, n);
      }
    }
    else
    {
      TreesIterator treesIt(nrOfVertices+1);
      uint8_t* currTree;
      while ((currTree = treesIt.next()))
      {
        for (uint j = 1; j <= nrOfVertices; j++)
        {
          automaton[cyclePos+j] = currTree[j] + cyclePos;
        }
        attachTreesToCycle(level+1, automaton, necklace, cycle, n);
      }
    }
  }

  void generatePartitions()
  {
    uint8_t aut[N]{};
    uint a[N]{};
    int k = 0;
    a[k] = N;
    bool generateTwoHalves = false;
    while(1)
    {
      if (isNumberOfCyclesInRange(k+1)){
        if (a[0] == N)
        {
          hasMoreThanHalfVs = true;
          automataSize = N;
          for (uint j=1; j < N; j++)
          {
            generatePartitions(N-j, j);
          }

          uint8_t lastAutomaton[N]{};
          for (uint k=0; k < N; k++)
            lastAutomaton[k] = (k+1)%N;
          reportAutomaton(lastAutomaton, N);
          autCnt++;
        }

        else
        {
          uint cnt[N]{};
          uint combInd[N]{};
          uint repInd[N]{};
          for (int i = 0; i <= k; i++)
            cnt[a[i]]++;
          for (uint i = 2; i <= N/2; i++)
          {
            if (cnt[i] > 1)
            {
              if (i == N/2)
                generateTwoHalves = true;
              else
              {
                uint64_t connAut = connautCnt[i-1];
                std::vector<uint> data(connAut);
                std::iota(std::begin(data), std::end(data), connAutUpToN[i-2]);
                std::vector<uint> res;
                uint chosen[N]{};
                chooseWithRep(chosen, 0, cnt[i], 0, data, res);
                combs[i] = std::move(res);
              }
            }
          }

          if (a[0] > N/2)
          {
            hasMoreThanHalfVs = true;
            automataSize = a[0];
            if (a[0] == N-1)
              automataTail[0] = N-1;
            else
              generateAutomata(0, k-1, 0, a+1, aut, combInd, repInd, cnt);

            if (N-a[0] == (uint)k)
            {
              for (uint j=1; j < automataSize; j++)
              {
                generatePartitions(automataSize-j, j);
              }

              uint8_t lastAutomaton[N]{};
              for (uint k=0; k < automataSize; k++)
                lastAutomaton[k] = (k+1)%automataSize;
              reportAutomatonWithTails(lastAutomaton);
              tailCnt = 0;
            }
          }

          else
          {
            hasMoreThanHalfVs = false;
            if (generateTwoHalves)
            {
              generateTwoHalves = false;
              uint start = connAutUpToN[N/2-2];
              uint end = connAutUpToN[N/2-1];

              for (uint i = start; i < end; i++)
              {
                for (uint j = i; j < end; j++)
                {
                  for (uint k = 0; k < 2; k++)
                  {
                    for (uint l = 0; l < N/2; l++)
                    {
                      if (k==0)
                        aut[l] = (*this)[i][l];
                      else
                        aut[N/2+l] = (*this)[j][l] + N/2;
                    }
                  }
                  if (N%2 != 0)
                    aut[N-1] = N-1;
                  reportAutomaton(aut, N);
                  autCnt++;
                }
              }
            }
            else
              generateAutomata(0, k, 0, a, aut, combInd, repInd, cnt);
          }
          for (uint i = 2; i <= N/2; i++)
            combs[i].clear();
        }
      }

      int remVal = 0;

      while (k >= 0 && a[k] == 1)
      {
        remVal += a[k];
        k--;
      }

      if (k < 0)  
        break;

      a[k]--;
      remVal++;

      while ((uint)remVal > a[k])
      {
        a[k+1] = a[k];
        remVal = remVal - a[k];
        k++;
      }

      a[k+1] = remVal;
      k++;
    }
    std::cout << "Number of automata generated: " << autCnt-1 << std::endl;
  }

  void generateAutomata(uint level, uint maxLevel, uint offset, uint* partition,
    uint8_t* automaton, uint* combInd, uint* repInd, uint* cnt)
  {
    if (level > maxLevel)
    {
      if (hasMoreThanHalfVs)
      {
        for (uint i = 0; i < N-automataSize; i++)
          automataTail[tailCnt*(N/2) + i] = automaton[i] + automataSize;
        tailCnt++;
      }
      else
      {
        reportAutomaton(automaton, N);
        autCnt++;
      }
      return;
    }

    uint n = partition[level];
    if (n == 1)
    {
      automaton[offset] = offset;
      generateAutomata(level+1, maxLevel, offset+1, partition, automaton, combInd, repInd, cnt);
    }
    else
    {
      if (!(combs[n].empty()))
      {
        const auto& nCombs = combs[n];
        uint nCnt = cnt[n];
        if (repInd[n] == 0)
        {
          for (uint i = 0; i < nCombs.size()/nCnt; i++)
          {
            for (uint j = 0; j < n; j++)
              automaton[offset+j] = (*this)[nCombs[nCnt*i]][j] + offset;
            repInd[n]++;
            generateAutomata(level+1, maxLevel, offset+n, partition, automaton, combInd, repInd, cnt);
            repInd[n]--;
            combInd[n]++;
          }
          combInd[n] = 0;
        }
        else
        {
          for (uint j = 0; j < n; j++)
            automaton[offset+j] = (*this)[nCombs[nCnt*combInd[n]+repInd[n]]][j] + offset;
          repInd[n]++;
          generateAutomata(level+1, maxLevel, offset+n, partition, automaton, combInd, repInd, cnt);
          repInd[n]--;
        }
      }
      else
      {
        for (uint64_t i = connAutUpToN[n-2]; i < connAutUpToN[n-1]; i++)
        {
          for (uint j = 0; j < n; j++)
            automaton[offset+j] = (*this)[i][j] + offset;
          generateAutomata(level+1, maxLevel, offset+n, partition, automaton, combInd, repInd, cnt);
        }
      }
    }
  }

  void chooseWithRep(uint* a, uint chosen, uint len, uint pos,
    const std::vector<uint>& data, std::vector<uint>& res)
  {
    if (chosen == len)
    {
      for (uint i = 0; i < len; i++)
        res.push_back(data[a[i]]);
      return;
    }
    for (uint i = pos; i < data.size(); i++)
    {
      a[chosen] = i;
      chooseWithRep(a, chosen + 1, len, i, data, res);
    }
  }

  bool isNumberOfCyclesInRange(uint cyclesNr)
  {
    return (cyclesNr <= maxNrOfCycles) && (cyclesNr >= minNrOfCycles);
  }

  void generate()
  {
    generatePartitions();
  }

  ~unaryAutomata() { delete[] automataTail; delete[] automata; }
};

#endif