/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  modified for WGD by Charles-Elie Rabier
  2011


  Newick tree reading/writing

=============================================================================*/

// c++ headers
#include <assert.h>
#include <stdio.h>

// spidir headers
#include "Tree.h"
#include "newick.h"
#include "parsing.h"
#include "logging.h"



namespace spidir {


float readDist(FILE *infile, int &depth)
{
    float dist = 0;
    fscanf(infile, "%f", &dist);
    return dist;
}


char readChar(FILE *stream, int &depth)
{
    char chr;
    do {
        if (fread(&chr, sizeof(char), 1, stream) != 1) {
            // indicate EOF
            return '\0';
        }
    } while (chr == ' ' || chr == '\n');
    
    // keep track of paren depth
    if (chr == '(') depth++;
    if (chr == ')') depth--;
    
    return chr;
}


char readUntil(FILE *stream, string &token, const char *stops, int &depth)
{
    char chr;
    token = "";
    while (true) {
        chr = readChar(stream, depth);
        if (!chr)
            return chr;
        
        // compare char to stop characters
        for (const char *i=stops; *i; i++) {
            if (chr == *i)
                return chr;
        }
        token += chr;
    }
}


int nodeNameCmp(const void *_a, const void *_b)
{
    Node *a = *((Node**) _a);
    Node *b = *((Node**) _b);
    
    if (a->nchildren == 0) {
        if (b->nchildren == 0)
            return 0;
        else
            return -1;
    } else {
        if (b->nchildren == 0)
            return 1;
        else
            return 0;
    }
}

Node *readNewickNode(FILE *infile, Tree *tree, Node *parent, int &depth)
{
  char  char1;
  Node *node;
  string token;
  bool wgdtest = false;
  // read first character
  if (!(char1  = readChar(infile, depth))) {
    printError("unexpected end of file");
    return NULL;
  }
    
    
  if (char1 == '(') {
    // read internal node
    
    int depth2 = depth;
    node = tree->addNode(new Node());
    if (parent)
      parent->addChild(node);
        
    // read all child nodes at this depth
    while (depth == depth2) {
      Node *child = readNewickNode(infile, tree, node, depth);
      if (!child)
	return NULL;
    }
        
    //read distance for this node
    char chr = readUntil(infile, token, "):,;", depth);

    if (chr==';'){
      return node;}
   
    if (!(token.length()<3)){ 
	  
      if (token.substr(0,3)=="WGD"){//found a WGD
	wgdtest=true;
      }
    }

    if (!wgdtest){ //no whole genome duplication
      node->longname = trim(token.c_str());
      node->dist = readDist(infile, depth);
      
      if (!(chr = readUntil(infile, token, "):,;", depth))){
	return NULL;}
      
    } 
    else{ //whole genome duplication
	  //we have to add a node
      Node *node2;
	 
      node2 = tree->addNode(new Node());	  
      node2->addChild(node);
	 
      // node was the child of parent created last
      parent->children[parent->nchildren - 1] = node2;
      node2->parent = parent;     
      node2->dist = readDist(infile, depth);
      node2->longname="WGD_before";      
      node->longname="WGD_at";
      node->dist=0;
      
      tree->settheWGD((tree->nWGD)+1);
	  
      //create a WGD to add to the tree
      WGDparam *WGDnew=new WGDparam(node2,node,
				    node->children[0],
				    atof(token.substr(3,token.length()-3).c_str()), 
				    node2->dist + node->children[0]->dist);
      
      tree->addWGD(WGDnew);
      
      if (!(chr = readUntil(infile, token, "):,", depth))){
	return NULL;}
    }

    return node;


  } else {
    // read leaf
            
    node = tree->addNode(new Node());
    if (parent)
      parent->addChild(node);
        
    char chr = readUntil(infile, token, "):,", depth);
        
    node->longname = char1 + trim(token.c_str());


    if (chr == ':') {
      node->dist = readDist(infile, depth);
      if (!(chr = readUntil(infile, token, ":),", depth)))
	return NULL;
    }

    return node;
  }
}


Tree *readNewickTree(FILE *infile, Tree *tree)
{
    bool newTree = false;
    if (!tree) {
        tree = new Tree();
        newTree = true;
    }

    int depth = 0;
    tree->root = readNewickNode(infile, tree, NULL, depth);
   

    // renumber nodes in a valid order
    // 1. leaves come first
    // 2. root is last

    if (tree->root != NULL) {
        tree->nodes[tree->root->name] = tree->nodes[tree->nnodes-1];    
        tree->nodes[tree->nnodes-1] = tree->root;
        tree->nodes[tree->root->name]->name = tree->root->name;
        tree->root->name = tree->nnodes - 1;
        
        qsort((void*) tree->nodes.get(), tree->nodes.size(), sizeof(Node*),
               nodeNameCmp);
        
        // update names
        for (int i=0; i<tree->nnodes; i++)
	  tree->nodes[i]->name = i;
        
        return tree;
    } else {
        if (newTree)
            delete tree;
        return NULL;
    }
}


Tree *readNewickTree(const char *filename, Tree *tree)
{
    bool newTree = false;
    if (!tree) {
        tree = new Tree();
        newTree = true;
    }

    FILE *infile = NULL;
    
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'\n", filename);
        if (newTree)
            delete tree;
        return NULL;
    }

    Tree *tree2 = readNewickTree(infile, tree);
    
    if (!tree2) {
        printError("format error in '%s'\n", filename);
        if (newTree)
            delete tree;
    }
    
    fclose(infile);
    
    return tree2;
}



// write out the newick notation of a tree
void writeNewickNode(FILE *out, Node *node, int depth, bool oneline)
{
    if (node->nchildren == 0) {
        if (!oneline)
            for (int i=0; i<depth; i++) fprintf(out, "  ");
        fprintf(out, "%s:%f", node->longname.c_str(), node->dist);
    } else {
        // indent
        if (oneline) {
            fprintf(out, "(");
        } else {
            for (int i=0; i<depth; i++) fprintf(out, "  ");
            fprintf(out, "(\n");
        }
        
        for (int i=0; i<node->nchildren - 1; i++) {
            writeNewickNode(out, node->children[i], depth+1, oneline);
            if (oneline)
                fprintf(out, ",");
            else
                fprintf(out, ",\n");
        }
        
        writeNewickNode(out, node->children[node->nchildren-1], 
                        depth+1, oneline);
        if (!oneline) {
            fprintf(out, "\n");            
            for (int i=0; i<depth; i++) fprintf(out, "  ");
        }
        fprintf(out, ")");
        
        if (depth > 0)
            fprintf(out, "%s:%f", node->longname.c_str(), node->dist);
        else
            fprintf(out, "%s", node->longname.c_str());
    }
}


// write out the newick notation of a tree
void writeNewickTree(FILE *out, Tree *tree, int depth, bool oneline)
{
    writeNewickNode(out, tree->root, 0, oneline);
    if (oneline)
        fprintf(out, ";");
    else
        fprintf(out, ";\n");
}


bool writeNewickTree(const char *filename, Tree *tree, bool oneline)
{
    FILE *out = NULL;
    
    if ((out = fopen(filename, "w")) == NULL) {
        printError("cannot write file '%s'\n", filename);
        return false;
    }

    writeNewickTree(out, tree, 0, oneline);
    fclose(out);
    return true;
}





} // namespace spidir
