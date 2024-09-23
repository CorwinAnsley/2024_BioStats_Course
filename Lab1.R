x = 1
y = 10
factor1 = 666
factor2 = 777

#print(factor1*factor2)

sum_to = function(target){
  i = 0
  sum = 0
  
  while (i < target+1) {
    sum = sum + i
    i = i + 1
  }
  return(sum)
}

a = sum_to(10)

b = 145000/5
c = (11 * 458) - (34667 / 444)

x + y

numbers1 = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
numbers2 = c(22, 24, 26, 28, 30)
numbers = c(numbers1, numbers2)
numbers_stats = c(mean(numbers), median(numbers), sd(numbers))

numbers[3]
index12 = numbers[12]
numbers[666]

numbers[5] = 999
numbers[6] = "potato"
numbers[7] = TRUE

vector_for_loop = c(101:200)

for (index in 1:100) 
{
  print(index)
}

for (index in 1:100) 
{
  print(vector_for_loop[index])
}

x = 1
if (x > 0) {print("I am bigger than zero")}
if (x < 0) {print("I am smaller than zero")}
if (x != 0) {print("I am not zero")}
if (x == 0) {print("I am zero")}

for (number in 1:100)
{
  if (number > 50)
  {
    print(number)
  }
}

print_only_even = function(vec) {
  for(num in vec) {
    if (num%%2 == 0) {
      print(num)
    }
  }
}

print_only_even(1:30)

