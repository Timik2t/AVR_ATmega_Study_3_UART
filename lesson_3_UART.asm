.include "m8515def.inc"

.def temp1=R17
.def temp2=R18
.def temp3=R19
.def counter1=r20
.def counter2=r21
.def counter2=r25
.def chislo1=r22
.def chislo2=r23
.def flag=r24

 rjmp main

.org $09          
 rjmp priem
 rjmp ozjidane

.org $0B

  in flag,sreg
  push flag 
  ldi temp1,$00       
 out DDRA,temp1
   clr temp1
   out DDRB,temp1                
 ldi temp2,$90      
 out UCSRB,temp2
 out sreg,flag
 pop flag

 reti

main:
bset 7
 ldi r30,low(RAMEND)
 out SPL,r30
 ldi r30,high(RAMEND)
 out SPH,r30

 clr temp1
 clr temp2
 clr counter1
 clr counter2
 ser chislo1
 ser chislo2
 clr flag
 clr temp3

 out DDRB,temp1
 ldi temp1,$FF        
 out DDRA,temp1


ldi temp3,$01
out ubrrh,temp3
ldi temp3,0
out ubrrl,temp3

 ldi temp1,$24       
 out UCSRC,temp1

 ldi temp2,$90      
 out UCSRB,temp2
 ser temp1
 out porta,temp1

 rjmp ozjidane

 .macro delay
	ser counter1
	delay1:
		ser counter2
		delay2:
			dec counter2
		brne delay2
		dec counter1
	brne delay1
.endmacro

 ozjidane:
   sbis pinb,7
  rjmp sravnenie1
   rjmp ozjidane
  
  sravnenie1:
  delay
  in temp1,pinb
  cpi temp1,$FF
   breq vvod
  rjmp sravnenie1

  vvod:
  sbi PINB,7
  in chislo1,PINB
   cpi chislo1,$FF
   brne sravnenie2
  rjmp vvod

  sravnenie2:
  delay
  in temp1,pinb
  cpi temp1,$FF
  breq otpravka
  rjmp sravnenie2

  otpravka:
  ser temp1
   out DDRB,temp1  
   ldi temp1,$48   
   out UCSRB,temp1
   out UDR,chislo1
   rjmp ozjidane1

ozjidane1: rjmp ozjidane1

  priem:
   in flag,sreg
   push flag  
   ldi temp1,$FF        
 out DDRA,temp1
   in chislo2,UDR
   out Porta,chislo2
   sbi Porta,7
   pop flag
   out sreg,flag
   delay
  reti


