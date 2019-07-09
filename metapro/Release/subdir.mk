################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../artp.cpp \
../lancaster.cpp \
../main.cpp \
../ordmeta.cpp \
../wFisher.cpp \
../wZ.cpp 

OBJS += \
./artp.o \
./lancaster.o \
./main.o \
./ordmeta.o \
./wFisher.o \
./wZ.o 

CPP_DEPS += \
./artp.d \
./lancaster.d \
./main.d \
./ordmeta.d \
./wFisher.d \
./wZ.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

main.o: ../main.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"main.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


