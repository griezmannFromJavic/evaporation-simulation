#include <stdio.h>

int main() {
    double num1, num2;
    char op;

    printf("Enter operator: \n");
    scanf("%c", &op);
    printf("Enter first number: \n");
    scanf("%lf", &num1);
    printf("Enter second number: \n");
    scanf("%lf", &num2);

    if (op == '+') {
        printf("%f", num1 + num2);
    } else if (op == '-') {
        printf("%f", num1 - num2);
    } else {
        printf("error");
    }

    return 0;
}
