#! /usr/bin/env python

if __name__ == '__main__':
    data = []
    with open("dataset.txt") as f:
        for i in f:
            data.append(i.split())
    data = map(float, data)
    print max(data), min(data)
