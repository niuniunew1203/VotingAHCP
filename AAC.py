import csv


def AAC(given_string):
    allowed_string = "AMVFIHYLWCDKSTGQEPRN"
    Dict = {char: 0 for char in allowed_string}

    for i in given_string:
        if i not in allowed_string:
            print("Invalid string")
            exit()
        else:
            Dict[i] += 1

    for i in Dict:
        Dict[i] = Dict[i]/len(given_string)

    with open('AAC.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header_row = ['Sequence'] + sorted(list(Dict.keys()))
        writer.writerow(header_row)
        data_row = [given_string] + [Dict[key] for key in sorted(Dict.keys())]
        writer.writerow(data_row)


given_string = "VSFAIKWEYVLLL"
AAC(given_string)
