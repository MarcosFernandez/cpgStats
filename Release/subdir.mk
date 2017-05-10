################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../common.c \
../counts.c \
../intersection.c \
../main.c \
../methBed.c \
../parsArgs.c \
../parseInput.c 

OBJS += \
./common.o \
./counts.o \
./intersection.o \
./main.o \
./methBed.o \
./parsArgs.o \
./parseInput.o 

C_DEPS += \
./common.d \
./counts.d \
./intersection.d \
./main.d \
./methBed.d \
./parsArgs.d \
./parseInput.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


