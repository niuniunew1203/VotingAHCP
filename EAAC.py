# def EAAC(given_string):

#     # _ should not be in the end csv file
#     allowed_string = "AMVFIHYLWCDKSTGQEPRN"
#     if (len(given_string) > 36):
#         print("Invalid string")
#         exit()
#     elif (len(given_string) < 36):
#         for i in given_string:
#             if i not in allowed_string:
#                 print("Invalid string")
#                 exit()
#         given_string = given_string.ljust(36, '_')
#         allowed_string += '_'
#         i = 0
#         j = i+5
#         while j <= len(given_string):
#             Dict = {char: 0 for char in allowed_string}
#             for k in given_string[i:j]:
#                 Dict[k] += 1

#             for m in Dict:
#                 Dict[m] = Dict[m]/5
#             # print(f'{i} Iteration :{given_string[i:j]}')
#             print(f'SW.{i+1}')
#             print(Dict)
#             i += 1
#             j += 1
#     else:
#         for i in given_string:
#             if i not in allowed_string:
#                 print("Invalid string")
#                 exit()
#         i = 0
#         j = i+5

#         while j <= len(given_string):
#             Dict = {char: 0 for char in allowed_string}
#             for k in given_string[i:j]:
#                 Dict[k] += 1

#             for m in Dict:
#                 Dict[m] = Dict[m]/5

#             # print(f'{i} Iteration: {given_string[i:j]}')
#             print(f'SW.{i+1}')
#             print(Dict)
#             i += 1
#             j += 1


# import csv


# def create_csv(given_string):
#     allowed_string = "AMVFIHYLWCDKSTGQEPRN"
#     if len(given_string) > 36:
#         print("Invalid string")
#         return
#     elif len(given_string) < 36:
#         for char in given_string:
#             if char not in allowed_string:
#                 print("Invalid string")
#                 return
#         given_string = given_string.ljust(36, ' ')
#         allowed_string = allowed_string.replace('_', '')
#     else:
#         for char in given_string:
#             if char not in allowed_string:
#                 print("Invalid string")
#                 return

#     rows = []
#     headers = ["Seq"]
#     i = 0
#     j = i + 5
#     while j <= len(given_string):
#         window_name = f"SW.{i+1}"
#         header = [f"{window_name}_{char}" for char in allowed_string]
#         headers.extend(header)

#         freq_dict = {char: 0 for char in allowed_string}
#         for k in given_string[i:j]:
#             if k != ' ':
#                 freq_dict[k] += 1

#         freq_list = [freq_dict[char]/5 for char in allowed_string]
#         row = [given_string.replace(' ', '')] + freq_list
#         rows.append(row)

#         i += 1
#         j += 1

#     filename = "EAAC.csv"
#     with open(filename, 'w', newline='') as csvfile:
#         writer = csv.writer(csvfile)
#         writer.writerow(headers)
#         writer.writerows(rows)

#     print(f"CSV file '{filename}' created successfully!")

# import csv


# def create_csv(given_string):
#     allowed_string = "AMVFIHYLWCDKSTGQEPRN"
#     if len(given_string) > 36:
#         print("Invalid string")
#         return
#     elif len(given_string) < 36:
#         for char in given_string:
#             if char not in allowed_string:
#                 print("Invalid string")
#                 return
#         given_string = given_string.ljust(36, ' ')
#         allowed_string = allowed_string.replace('_', '')
#     else:
#         for char in given_string:
#             if char not in allowed_string:
#                 print("Invalid string")
#                 return

#     headers = ["Seq"]
#     row = [given_string.replace(' ', '')]

#     i = 0
#     j = i + 5
#     while j <= len(given_string):
#         window_name = f"SW.{i+1}"
#         header = [f"{window_name}_{char}" for char in allowed_string]
#         headers.extend(header)

#         freq_dict = {char: 0 for char in allowed_string}
#         for k in given_string[i:j]:
#             if k != ' ':
#                 freq_dict[k] += 1

#         freq_list = [freq_dict[char]/5 for char in allowed_string]
#         row.extend(freq_list)

#         i += 1
#         j += 1

#     data = [row]

#     filename = "EAAC.csv"
#     with open(filename, 'w', newline='') as csvfile:
#         writer = csv.writer(csvfile)
#         writer.writerow(headers)
#         writer.writerows(data)

#     print(f"CSV file '{filename}' created successfully!")

import csv


def EAAC(given_string):
    allowed_string = "AMVFIHYLWCDKSTGQEPRN"
    if len(given_string) > 36:
        print("Invalid string")
        return
    elif len(given_string) < 36:
        for char in given_string:
            if char not in allowed_string:
                print("Invalid string")
                return
        given_string = given_string.ljust(36, ' ')
        allowed_string = allowed_string.replace('_', '')
    else:
        for char in given_string:
            if char not in allowed_string:
                print("Invalid string")
                return

    headers = ["Seq"]
    row = [given_string.replace(' ', '')]

    i = 0
    j = i + 5
    while j <= len(given_string):
        window_name = f"SW.{i+1}"
        freq_dict = {char: 0 for char in allowed_string}
        for k in given_string[i:j]:
            if k != ' ':
                freq_dict[k] += 1

        freq_list = [freq_dict[char]/5 for char in sorted(allowed_string)]
        headers.extend(
            [f"{window_name}_{char}" for char in sorted(allowed_string)])
        row.extend(freq_list)

        i += 1
        j += 1

    data = [row]

    filename = "EAAC.csv"
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(data)

    print(f"CSV file '{filename}' created successfully!")


seq = 'VSFAIKWEYVLLL'
EAAC(seq)
