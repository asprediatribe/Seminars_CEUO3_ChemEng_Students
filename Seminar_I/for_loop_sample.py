# this script describes the main components of a for loop in Python:
# "range ()": Python's range() function returns an
# iterable list of values starting at zero and ending at n-1. For example, when range(3) is called, the values 0, 1,
# 2 are returned. 3 is not part of the output, even though the function input was range(3).
# "len()" : The Python’s  len() function  returns the length of an object. For example, it can return the number of items in a list
#You can use the function with many different data types.

ion_list = ['Na', '+1', 'cation']
for i in ion_list:
    print(i)

# ion_list = ['Na', '+1', 'cation']
# for i in range(3):
#     print(ion_list[i])

#ion_list = ['Na', '+1', 'cation']
#for i in range(len(ion_list)):
#    print(ion_list[i])
