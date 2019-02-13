#ifndef UNARYAUTOMATA_TREES_HPP_
#define UNARYAUTOMATA_TREES_HPP_

#include "constants.hpp"

struct TrieNode
{
  uint n = 0;
  TrieNode** children;

  TrieNode(uint nextSize)
  {
    n = nextSize;
    children = new TrieNode*[nextSize];
    
    for (uint i = 0; i < n; i++)
      children[i] = nullptr;
  }

  TrieNode()
  {
    children = nullptr;
  }

  TrieNode(const TrieNode&) = delete;
  TrieNode& operator=(const TrieNode&) = delete;

  ~TrieNode()
  {
    if (children)
    {
      for (uint i = 0; i < n; i++)
      {
        if (children[i])
          delete children[i];
      }
      delete[] children;
    }
  }
};

void insertConf(uint* conf, uint* necklace, uint n, TrieNode* trie)
{
  for (uint i = 0; i < n; i++)
  {
    uint ind = conf[i] - treesNrUpToN[necklace[i]];
    if (!trie->children[ind])
    {
      if (i == n-1)
        trie->children[ind] = new TrieNode();
      else
        trie->children[ind] = new TrieNode(treesNr[necklace[i+1]+1]);
    }
    trie = trie->children[ind];
  }
}

bool isConfInTrie(uint* conf, uint* necklace, uint n, TrieNode* trie)
{
  for (uint i = 0; i < n; i++)
  {
    uint ind = conf[i] - treesNrUpToN[necklace[i]];
    if (!trie->children[ind])
      return false;
    else
      trie = trie->children[ind];
  }
  return true;
}

struct TreesIterator
{
  uint N;
  uint* predIndex;
  uint* lastPred;
  uint8_t* currTree;
  uint pos;
  bool firstTree;

  TreesIterator(uint n)
  : N(n),
  predIndex(new uint[n]),
  lastPred(new uint[n]),
  currTree(new uint8_t[n]),
  pos(n-1),
  firstTree(true)
  {}

  uint8_t* next()
  {
    if (firstTree)
    {
      firstTree = false;
      currTree[0] = 0;
      predIndex[0] = 1;
      lastPred[0] = 0;
      for (uint i = 1; i < N; i++)
      {
        currTree[i] = i-1;
        predIndex[i] = i+1;
        lastPred[i] = 0;
      }
      return currTree;
    }

    if (pos <= 1)
      return nullptr;

    currTree[pos]--;
    if (pos < N-1 && (currTree[pos] != 0 || currTree[pos-1] != 0))
    {
      uint diff = pos - predIndex[currTree[pos]];
      while (pos < N-1)
      {
        lastPred[pos] = predIndex[currTree[pos]];
        predIndex[currTree[pos]] = pos;
        pos++;
        currTree[pos] = currTree[pos - diff];
      }
    }
    while (currTree[pos] == 0 && pos > 0)
    {
      pos--;
      predIndex[currTree[pos]] = lastPred[pos];
    }
    return currTree;
  }

  ~TreesIterator()
  {
    delete[] currTree;
    delete[] lastPred;
    delete[] predIndex;
  }
};

template <uint N>
class TreesUpToN
{
public:
  uint predIndex[N];
  uint lastPred[N];
  uint8_t* trees;

  TreesUpToN() : trees(new uint8_t[treesNrUpToN[N] * N])
  {
    for (uint i = 1; i <= N; i++)
    {
      generateTreesWithNVertices(i, treesNrUpToN[i-1]);
    }
  }

  inline uint8_t* operator[](const uint64_t k) { return trees + N * k; }
  inline const uint8_t* operator[](const uint64_t k) const { return trees + N * k; }

  void printTree(uint64_t nr)
  {
    for (uint i = 0; i < N; i++)
    {
      std::cout << +(*this)[nr][i] << " ";
    }
  }

  void printTrees()
  {
    for (uint64_t i = 0; i < treesNrUpToN[N]; i++)
    {
      std::cout << i+1 << ": ";
      printTree(i);
      std::cout << std::endl;
    }
  }

  void generateTreesWithNVertices(const uint n, const uint64_t treeNr)
  {
    trees[treeNr*N] = 0;
    predIndex[0] = 1;
    lastPred[0] = 0;
    for (uint i = 1; i < n; i++)
    {
      trees[treeNr*N + i] = i-1;
      predIndex[i] = i+1;
      lastPred[i] = 0;
    }
    for (uint i = n; i < N; i++)
    {
      trees[treeNr*N + i] = 0;
    }

    uint8_t* currTree = (*this)[treeNr+1];
    uint pos = n-1;
    while (pos > 1)
    { 
      std::memcpy(currTree, currTree - N, N*sizeof(uint8_t));
      currTree[pos]--;
      if (pos < n-1 && (currTree[pos] != 0 || currTree[pos-1] != 0))
      {
        uint diff = pos - predIndex[currTree[pos]];
        while (pos < n-1)
        {
          lastPred[pos] = predIndex[currTree[pos]];
          predIndex[currTree[pos]] = pos;
          pos++;
          currTree[pos] = currTree[pos - diff];
        }
      }
      while (currTree[pos] == 0 && pos > 0)
      {
        pos--;
        predIndex[currTree[pos]] = lastPred[pos];
      }
      currTree += N;
    }
  }

  ~TreesUpToN() { delete[] trees; }
};

#endif