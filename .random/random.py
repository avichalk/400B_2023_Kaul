s = 'violet,cyan,blue,green,yellow,orange,red'
arr = s.split(',')

string = 'Centeralized Supermassive Black Holes!'
out = ''

x = 0
for i in string:
    if x == len(arr)-1:
        x = 0
    out += '{\color{' + arr[x] + '}' + i + '}'
    x += 1

print(out)
