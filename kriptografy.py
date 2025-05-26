from epileptic import *
from optimus import *
from Crypto.Cipher import AES
from Crypto.Hash import SHA256

def str_to_bytes(s):
    return s.encode('utf-8')

def bytes_to_str(b):
    return b.decode('utf-8')

def bytes_to_int(b):
    return int.from_bytes(b, byteorder='big')

def int_to_bytes(i):
    return i.to_bytes((i.bit_length() + 7) // 8, byteorder='big')


class PKI:
    def __init__(self):
        self.users = {}


    def register_user(self, name, pk, alg):
        if name in self.users:
            self.users[name][alg] = pk
        else:
            self.users[name] = {alg: pk}
 

    def get_public_key(self, name, alg):
        return self.users[name][alg]
    



class BaseUser:
    def __init__(self, name, pki):
        self.EC = P224()
        self.secret_key, self.public_key = self._ec_p224_keygen()
        self.name = name
        self.pki = pki

    def _ec_p224_keygen(self):
        sk = rand_less_than(self.EC.n - 1)
        pk = sk * self.EC.G

        return sk, pk



class UserSecretSharer(BaseUser):
    def __init__(self, name, pki):
        super().__init__(name, pki)
        self.shared_secrets = {}
        self.pki.register_user(name, self.public_key, "DH")


    def send_DH_init(self, name):
        other_pk = self.pki.get_public_key(name)
        shared_secret = SHA256.new((self.secret_key * other_pk).to_bytes()).digest()
        self.shared_secrets[name] = shared_secret

        return {
            "package_type": "DH_init",
            "name": self.name,
        }


    def recieve_DH_init(self, request):
        other_pk = self.pki.get_public_key(request["name"])
        shared_secret = SHA256.new((self.secret_key * other_pk).to_bytes()).digest()
        self.shared_secrets[request["name"]] = shared_secret    

    


class UserHybridTalker(BaseUser):
    def __init__(self, name, pki):
        super().__init__(name, pki)
        self.chats = {}
        self.pki.register_user(name, self.public_key, "DHE")


    def hybrid_send(self, name, message):
        if name in self.chats:
            k = self.chats[name]["k"]
            Cm = AES.new(k, AES.MODE_CBC).encrypt(str_to_bytes(message))
            self.chats[name]["messages"].append(f"Sent:  {message}")

            return {
                "package_type": "established_send",
                "name": self.name,
                "message": Cm
            }
        

        else:
            k = int_to_bytes(rand_int(128))
            Cm = AES.new(k, AES.MODE_CBC).encrypt(str_to_bytes(message))
            other_pk = self.pki.get_public_key(name)
            esk, epk = self._ec_p224_keygen()
            encaps_key = SHA256.new((esk * other_pk).to_bytes()).digest()[0:16]
            Ck = AES.new(encaps_key, AES.MODE_CBC).encrypt(k)
            self.chats[name] = {
                "k": k,
                "messages": [f"Sent:  {message}"]
            }

            return {
                "package_type": "new_send",
                "name": self.name,
                "DHE_piece": epk.to_bytes(),
                "encapsulated_key": Ck,
                "message": Cm
            }


    def hybrid_receive(self, request):
        request_type = request["package_type"]

        if request_type == "new_send":
            other_pk = EllipticCurvePoint.from_bytes(request["DHE_piece"], self.EC)
            encaps_key = SHA256.new((self.secret_key * other_pk).to_bytes()).digest()[0:16]
            k = AES.new(encaps_key, AES.MODE_CBC).decrypt(request["encapsulated_key"])
            message = bytes_to_str(AES.new(k, AES.MODE_CBC).decrypt(request["message"]))

            self.chats[request["name"]] = {
                "k": k,
                "messages": [f"Received:  {message}"]
            }
        

        elif request_type == "established_send":
            k = self.chats[request["name"]]["k"]
            message = bytes_to_str(AES.new(k, AES.MODE_CBC).decrypt(request["message"]))
            self.chats[request["name"]]["messages"].append(f"Received:  {message}")



class DigitalSignatureUser(BaseUser):
    def __init__(self, name, pki):
        super().__init__(name, pki)
        self.pki.register_user(name, self.public_key, "ECDSA")


    def sign_message(self, message):
        z = SHA256.new(str_to_bytes(message)).digest()[0:28]

        s = 0
        while s == 0:
            r = 0
            while r == 0:
                k = rand_less_than(self.EC.n - 1)
                R = k * self.EC.G
                r = R.x % self.EC.n

            s = (pow(k, -1, self.EC.n) * (z + self.secret_key * r)) % self.EC.n
        
        signature = (r, s)

        return {
            "package_type": "signed_message",
            "name": self.name,
            "message": message,
            "signature": signature
        }
    

    def verify_signature(self, sign_package):
        pk = self.pki.get_public_key(sign_package["name"], "ECDSA")
        r, s = sign_package["signature"]
        z = SHA256.new(str_to_bytes(sign_package["message"])).digest()[0:28]
        
        if not ((0 < r < self.EC.n) and (0 < s < self.EC.n)):
            return False
        
        w = pow(s, -1, self.EC.n)
        u1 = (z * w) % self.EC.n
        u2 = (r * w) % self.EC.n

        R = (u1 * self.EC.G) + (u2 * pk)

        if R == self.EC.inf():
            return False
        
        return R.x % self.EC.n == r

