cat test.txt | java -jar CPSM.jar -p 2,3 -m 1 -b 4 --hw 20
#0       1       2       3       4       5
#1.8328  1.0465  0.9095  0.9969  1.8789  0

cat test.txt | java -jar CPSM.jar -p 2,3 -m 1 -b 4 --hw 20 --p2 3
#0       1       2       3       4       5
#2.4031  1.657   1.3167  1.436   2.8772  0
#4.0485  2.4372  2.067   2.2651  6.0804  0
#3       1.8625  1.423   1.3646  2.5787  0
