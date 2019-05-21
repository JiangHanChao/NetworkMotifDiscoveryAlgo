#!/usr/bin/env python3
# -*- coding:utf-8 -*-


# if array doesn't contain the node, it will insert the node into the array and return false, else return true
def Binarysort(array, node):
    low = 0
    high = len(array) - 1
    while low <= high:
        mid = int((low + high)/2)
        if node < array[mid]:
            high = mid - 1
        elif node > array[mid]:
            low = mid + 1
        else:
            return True
    array.append(node)
    return False


def Binarycontain(array, node):
    low = 0
    high = len(array) - 1
    while low <= high:
        mid = int((low + high)/2)
        if node < array[mid]:
            high = mid - 1
        elif node > array[mid]:
            low = mid + 1
        else:
            return mid
    return -1


def Getlabel(mat, nodelist, Subg):
    label = ""
    nodetype = ""
    for m in Subg:
        for n in Subg:
            label += str(mat[m][n])
        nodetype += str(nodelist[m].nodetype)

    label += "_" + nodetype
    return label


def Swap(array, i, j):
    if len(array) < 2:
        raise ValueError("swap array size < 1")
    if i == j:
        return
    t = array[i]
    array[i] = array[j]
    array[j] = t
