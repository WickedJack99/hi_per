import math

def calc():
    with open("output.txt", "w") as f:
        values = ""
        for i in range(0,256):
            value = fib(i)
            val = (math.sin(value) * math.tan(value) / math.sqrt(math.cos(value) + 2))
            values += f"{val}, "
        f.write(values)

def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

calc()