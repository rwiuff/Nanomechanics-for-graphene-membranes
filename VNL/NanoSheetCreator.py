sheet = NanoSheet(2,2,Carbon)

nlprint(sheet)
sheet = sheet.repeat(1,5,20)
sheet = sheet.center()
nlprint(sheet)

#nlsave('test',sheet)