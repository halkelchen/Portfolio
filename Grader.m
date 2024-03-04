clc, clear
for student = 1 : 10
exam1=input('Enter score for exam 1');
exam2=input('Enter score for exam 2');
exam3=input('Enter score for exam 3');
exam4=input('Enter score for exam 4');
assignments1=input('Enter score for assignment 1');
assignments2=input('Enter score for assignment 2');
assignments3=input('Enter score for assignment 3');
assignments4=input('Enter score for assignment 4');
A = [exam1 assignments1; exam2 assignments2; exam3 assignments3; exam4 assignments4];
Sum = sum(sum(A));
%need to fix above line
if A>=450 & A<=500
    grade = 'A';
elseif A>=351 & A<=449
  grade ='B';
elseif A>=301 & A<=350
  grade ='C';
elseif A==300
  grade ='D';
else
      grade ="F";
end
end
classavg = mean(A);
classstd = std(A);
if (grade >= classavg + 2*classstd)% I can't figure out how to use teh strcmp, it keeps giving me a comparison between doubles and strings error
    grade = 'A';
elseif  (grade >= classavg + 1*classstd) && (grade<= classavg + 2*classstd)
  grade ='B';
elseif (grade >= classavg && grade<= classavg + 1*classstd)
  grade ='C';
elseif (grade >= classavg - 1*classstd && grade<= classavg)
  grade ='D';
else
      grade ="F";
end






