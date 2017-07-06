#!/usr/bin/env python3

def main():
    with open('D:/Project/Rice/mature.fa') as mature, open('D:/Project/Rice/osa_mature.fa', 'w') as osa_mature:
        line_count = 0
        for line in mature:
            if line.startswith('>osa'):
                osa_mature.write(line + next(mature).replace('U', 'T'))
                line_count += 1
        print(line_count)

if __name__ == '__main__':
    main()
