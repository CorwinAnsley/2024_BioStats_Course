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
   