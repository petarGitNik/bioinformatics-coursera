# Regarding code-graded problems

To redirect a file to STDIN, you can use the '<' symbol. A basic example Python3 script can be found below:

```python3
import sys # you must import "sys" to read from STDIN
lines = sys.stdin.read().splitlines() # read in the input from STDIN
for i in xrange(len(lines)):
    print('Line ' + str(i+1) + ' is ' + str(len(lines[i])) + ' characters long.')
```

To run the script, enter the following command:

```
example_read_stdin.py < example_input.txt
```
