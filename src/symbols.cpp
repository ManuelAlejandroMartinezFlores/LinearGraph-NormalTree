#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>

using namespace std;

const double EPSILON_SYM = 1e-15;

bool is_zero(double x) {
    return abs(x) < EPSILON_SYM;
}

class Expression {
    public: 
    virtual ~Expression() = default;
    virtual double evaluate(unordered_map<string, double> context) = 0;
    virtual string to_str() = 0;
    virtual double to_abs() = 0;
    virtual bool is_scalar() = 0;
    virtual bool is_atom() = 0;
    virtual bool is_reciprocal() = 0;
    virtual bool is_product() = 0;
    virtual string to_str_no_coeff() {return to_str();}
    virtual bool isZero() {
        return false;
    }
};

struct CompareSharedPtrExpression {
    bool operator()(const shared_ptr<Expression>& lhs, shared_ptr<Expression>& rhs) {
        return lhs->to_str_no_coeff() > rhs->to_str_no_coeff();
    }
};


class Scalar: public Expression {
    public:
    double value;
    Scalar(double v) {value = v;}
    double evaluate(unordered_map<string, double> context) override {
        return value;
    }
    string to_str() override {
        return to_string(value);
    }
    double to_abs() override {return abs(value);}
    bool is_scalar() override {return true;}
    bool is_atom() override {return false;}
    bool is_product() override {return false;}
    bool is_reciprocal() override {return false;}
    bool isZero() override {return is_zero(value);}
};

class Symbol: public Expression {
    public:
    string symbol; 
    Symbol(string s) {symbol = s;}
    double evaluate(unordered_map<string, double> context) override {
        return context[symbol];
    }
    string to_str() override {
        return symbol;
    }
    double to_abs() override {return 99;}
    bool is_scalar() override {return false;}
    bool is_atom() override {return true;}
    bool is_product() override {return false;}
    bool is_reciprocal() override {return false;}
};


class Reciprocal: public Expression {
    public:
        shared_ptr<Expression> content;
        Reciprocal(shared_ptr<Expression> c) : content(c) {}
        double evaluate(unordered_map<string, double> context) override {
            return 1 / content->evaluate(context);
        }
        string to_str() override {
            return "(1/" + content->to_str() + ")";
        }
        double to_abs() override {return 99;}
        bool is_scalar() override {return false;}
        bool is_atom() override {return true;}
        bool is_product() override {return false;}
        bool is_reciprocal() override {return true;}
};


class Multiplication : public Expression {
    public:
        double scalar;
        vector<shared_ptr<Expression>> factors;
        Multiplication() : scalar(1) {}
        double evaluate(unordered_map<string, double> context) override {
            double ans = scalar;
            for (auto factor : factors) {
                ans *= factor->evaluate(context);
            }
            return ans;
        }
        string to_str() override {
            string ans = to_string(scalar) ;
            for (auto factor : factors) {
                ans += "*" + factor->to_str();
            }
            return ans;
        }
        double to_abs() override {
            return 99;
        }
        void add_factor(shared_ptr<Expression> factor) {
            factors.push_back(factor);
            push_heap(factors.begin(), factors.end(), CompareSharedPtrExpression());
        }
        bool is_scalar() override {return false;}
        bool is_atom() override {return false;}
        bool is_product() override {return true;}
        bool is_reciprocal() override {return false;}
        string to_str_no_coeff() override {
            string ans = "";
            if (factors.size() > 0) {
                ans += factors[0]->to_str();
            }
            for (int i=1; i<factors.size(); i++) {
                ans += "*" + factors[i]->to_str() ;
            }
            return ans;
        }
        bool isZero() override {
            return is_zero(scalar);
        }
};

class Sum : public Expression {
    public:
        double independent;
        vector<shared_ptr<Expression>> sumands;
        Sum() : independent(0) {}
        double evaluate(unordered_map<string, double> context) override {
            double ans = 0;
            for (auto sumand : sumands) {
                ans += sumand->evaluate(context);
            }
            return ans;
        }
        string to_str() override {
            string ans = "(" + to_string(independent) ;
            for (auto sumand : sumands) {
                ans += "+" + sumand->to_str();
            }
            return ans + ")";
        }
        double to_abs() override {
            return 99;
        }
        void add_sumand(shared_ptr<Expression> other) {
            if (other->is_scalar()) {
                auto sumand_scalar = dynamic_pointer_cast<Scalar>(other);
                independent += sumand_scalar->value;
            } else {
                // Simplify similar terms
                bool flag = false;
                vector<shared_ptr<Expression>> new_sumands;
                for (int i=0; i<sumands.size(); i++) {
                    auto sumand = sumands[i];
                    if (sumand->to_str_no_coeff() == other->to_str_no_coeff()) {
                        flag = true;
                        if (other->is_product()) {
                            // New sumand is product
                            auto other_prod = dynamic_pointer_cast<Multiplication>(other);
                            if (sumand->is_product()) {
                                // Similar term is product
                                auto temp = make_shared<Multiplication>();
                                auto sumand_prod = dynamic_pointer_cast<Multiplication>(sumand);
                                for (auto factor : other_prod->factors) {
                                    temp->add_factor(factor);
                                }
                                temp->scalar = other_prod->scalar + sumand_prod->scalar;
                                if (!is_zero(temp->scalar)) {
                                    new_sumands.push_back(temp);
                                }
                            } else if (sumand->is_atom()) {
                                // Similar term is atom
                                cout << "----" << endl;
                                auto temp = make_shared<Multiplication>();
                                for (auto factor : other_prod->factors) {
                                    temp->add_factor(factor); 
                                }
                                temp->scalar = other_prod->scalar + 1;
                                if (!is_zero(temp->scalar)) {
                                    new_sumands.push_back(temp);
                                }
                            }
                        } else if (other->is_atom()) {
                            // New sumand is atom 
                            if (sumand->is_product()) {
                                // Similar term is product
                                auto sumand_prod = dynamic_pointer_cast<Multiplication>(sumand);
                                auto temp = make_shared<Multiplication>();
                                for (auto factor : sumand_prod->factors) {
                                    temp->add_factor(factor); 
                                }
                                temp->scalar = sumand_prod->scalar + 1;
                                if (!is_zero(temp->scalar)) {
                                    new_sumands.push_back(temp);
                                }
                            } else if (sumand->is_atom()) {
                                // Similar term is atom 
                                auto temp = make_shared<Multiplication>();
                                temp->scalar = 2;
                                temp->add_factor(sumand);
                                new_sumands.push_back(temp);
                            }
                        }
                    } else {
                        new_sumands.push_back(sumand);
                    }
                }
                if (!flag) {
                    new_sumands.push_back(other);
                } 
                sumands = new_sumands;
                make_heap(sumands.begin(), sumands.end(), CompareSharedPtrExpression());
            }
        }
        bool is_scalar() override {return false;}
        bool is_atom() override {return false;}
        bool is_product() override {return false;}
        bool is_reciprocal() override {return false;}
        bool isZero() override {
            if (!is_zero(independent)) {
                return false;
            } for (auto sumand : sumands) {
                if (!sumand->isZero()) {
                    return false;
                }
            }
            return true;
        }
};




shared_ptr<Expression> operator*(shared_ptr<Expression> lhs, shared_ptr<Expression> other) {
    if (lhs->isZero()) {
        return make_shared<Scalar>(0);
    }
    if (lhs->is_scalar()) {
        // cout << "scalar" << endl;
        // Scalar mulitplication
        auto lhs_scalar = dynamic_pointer_cast<Scalar>(lhs);

        if (is_zero(lhs_scalar->value-1)) {
            return other;
        } 
        if (other->is_scalar()) {
            // With scalar
            auto other_scalar = dynamic_pointer_cast<Scalar>(other);
            return make_shared<Scalar>(lhs_scalar->value * other_scalar->value);
        } if (other->is_atom()) {
            // With symbol or reciprocal
            shared_ptr<Multiplication> ans = make_shared<Multiplication>();
            ans->scalar = lhs_scalar->value;
            ans->add_factor(other);
            return ans;
        } if (other->is_product()) {
            // With product
            auto other_product = dynamic_pointer_cast<Multiplication>(other);
            auto ans = make_shared<Multiplication>();
            ans->scalar = lhs_scalar->value * other_product->scalar;
            for (auto factor : other_product->factors) {
                ans->add_factor(factor);
            }
            return ans;
        } else {
            // With sum
            auto other_sum = dynamic_pointer_cast<Sum>(other);
            auto ans = make_shared<Sum>();
            ans->independent = other_sum->independent * lhs_scalar->value;
            for (auto sumand : other_sum->sumands) {
                ans->add_sumand(lhs_scalar * sumand);
            }
            if (is_zero(ans->independent) && ans->sumands.size()==1) {
                return ans->sumands[0];
            } if (ans->sumands.size() == 0) {
                return make_shared<Scalar>(ans->independent);
            }
            return ans;
        }
    } if (lhs->is_atom() && !lhs->is_reciprocal()) {
        // cout << "symbol" << endl;
        // Symbol multiplication
        auto lhs_symbol = dynamic_pointer_cast<Symbol>(lhs);
        if (other->is_scalar()) {
            // With scalar
            return other * lhs;
        } if (other->is_atom()) {
            if (other->is_reciprocal()) {
                // With reciprocal
                auto other_rec = dynamic_pointer_cast<Reciprocal>(other);
                if (lhs_symbol->symbol == other_rec->content->to_str()) {
                    return make_shared<Scalar>(1);
                }
            }
            // With symbol
            auto ans = make_shared<Multiplication>();
            ans->add_factor(other);
            ans->add_factor(lhs);
            return ans;
        } if (other->is_product()) {
            // With product
            auto other_product = dynamic_pointer_cast<Multiplication>(other);
            auto ans = make_shared<Multiplication>();
            ans->scalar = other_product->scalar;
            bool cancels = false;
            // Check if it cancels with a reciprocal
            for (auto factor : other_product->factors) {
                if (!cancels && factor->is_reciprocal()) {
                    auto factor_rec = dynamic_pointer_cast<Reciprocal>(factor);
                    if (lhs->to_str() == factor_rec->content->to_str()) {
                        cancels = true;
                        continue;
                    }
                }
                ans->add_factor(factor);
            }
            if (!cancels) {
                ans->add_factor(lhs);
            }
            return ans;
        } else {
            // With sum
            auto other_sum = dynamic_pointer_cast<Sum>(other);
            auto ans = make_shared<Sum>();
            for (auto sumand : other_sum->sumands) {
                ans->add_sumand(lhs * sumand);
            }
            ans->add_sumand(make_shared<Scalar>(other_sum->independent) * lhs);
            if (is_zero(ans->independent) && ans->sumands.size()==1) {
                return ans->sumands[0];
            } if (ans->sumands.size() == 0) {
                return make_shared<Scalar>(ans->independent);
            }
            return ans;
        }
    } if (lhs->is_reciprocal()) {
        // cout << "rec" << endl;
        // Reciprocal multiplication
        auto lhs_rec = dynamic_pointer_cast<Reciprocal>(lhs);
        if (other->is_scalar()) {
            // With scalar
            return other * lhs_rec;
        } if (other->is_atom()) {
            if (other->is_reciprocal()) {
                // With reciprocal
                auto ans = make_shared<Multiplication>();
                ans->add_factor(other);
                ans->add_factor(lhs);
                return ans;
            }
            // With symbol
            return other * lhs;
        } if (other->is_product()) {
            // With product
            auto other_product = dynamic_pointer_cast<Multiplication>(other);
            auto ans = make_shared<Multiplication>();
            ans->scalar = other_product->scalar;
            bool cancels = false;
            // Check if it cancels with a symbol
            for (auto factor : other_product->factors) {
                if (!cancels && !factor->is_reciprocal() && factor->is_atom()) {
                    auto factor_rec = dynamic_pointer_cast<Symbol>(factor);
                    if (lhs_rec->content->to_str() == factor_rec->symbol) {
                        cancels = true;
                        continue;
                    }
                }
                ans->add_factor(factor);
            }
            if (!cancels) {
                ans->add_factor(lhs);
            }
            return ans;
        } else {
            // With sum
            auto other_sum = dynamic_pointer_cast<Sum>(other);
            auto ans = make_shared<Sum>();
            for (auto sumand : other_sum->sumands) {
                ans->add_sumand(lhs_rec * sumand);
            }
            ans->add_sumand(make_shared<Scalar>(other_sum->independent) * lhs);
            if (is_zero(ans->independent) && ans->sumands.size()==1) {
                return ans->sumands[0];
            } if (ans->sumands.size() == 0) {
                return make_shared<Scalar>(ans->independent);
            }
            return ans;
        }
    } if (lhs->is_product()) {
        // cout << "product" << endl;
        // Product multiplication
        auto lhs_prod = dynamic_pointer_cast<Multiplication>(lhs);
        if (other->is_scalar()) {
            // With scalar
            return other * lhs_prod;
        } if (other->is_atom()) {
            // With symbol or reciprocal
            return other * lhs_prod;
        } if (other->is_product()) {
            // With product
            auto other_product = dynamic_pointer_cast<Multiplication>(other);
            auto ans = make_shared<Multiplication>();
            ans->scalar = other_product->scalar * lhs_prod->scalar;
            for (auto factor : other_product->factors) {
                ans->factors.emplace_back(factor);
            }
            shared_ptr<Expression> newans = ans;
            for (auto factor : lhs_prod->factors) {
                newans = (factor * newans);
            }
            return newans;
        } else {
            // With sum
            auto other_sum = dynamic_pointer_cast<Sum>(other);
            auto ans = make_shared<Sum>();
            for (auto sumand : other_sum->sumands) {
                ans->add_sumand(lhs_prod * sumand);
            }
            ans->add_sumand(make_shared<Scalar>(other_sum->independent) * lhs);
            if (is_zero(ans->independent) && ans->sumands.size()==1) {
                return ans->sumands[0];
            } if (ans->sumands.size() == 0) {
                return make_shared<Scalar>(ans->independent);
            }
            return ans;
        }
    } else {
        // Sum multiplication
        // cout << "sum" << endl;
        auto lhs_sum = dynamic_pointer_cast<Sum>(lhs);

        if (other->is_scalar()) {
            // With scalar
            return other * lhs_sum;
        } if (other->is_atom()) {
            // With symbol or reciprocal
            return other * lhs_sum;
        } if (other->is_product()) {
            // With product
            return other * lhs_sum;
        } else {
            // With sum
            auto other_sum = dynamic_pointer_cast<Sum>(other);
            auto ans = make_shared<Sum>();
            ans->independent = lhs_sum->independent * other_sum->independent;
            for (auto sumand1 : lhs_sum->sumands) {
                for (auto sumand2 : other_sum->sumands) {
                    ans->add_sumand(sumand1 * sumand2);
                }
                ans->add_sumand(make_shared<Scalar>(other_sum->independent) * sumand1);
            }
            for (auto sumand2 : other_sum->sumands) {
                ans->add_sumand(make_shared<Scalar>(lhs_sum->independent) * sumand2);
            }
            if (is_zero(ans->independent) && ans->sumands.size()==1) {
                return ans->sumands[0];
            } if (ans->sumands.size() == 0) {
                return make_shared<Scalar>(ans->independent);
            }
            return ans;
        }
    }
    return make_shared<Scalar>(-99);
}

shared_ptr<Expression> operator*(double lhs, shared_ptr<Expression> other) {
    return make_shared<Scalar>(lhs) * other;
}

shared_ptr<Expression> operator*(shared_ptr<Expression> lhs, double other) {
    return make_shared<Scalar>(other) * lhs;
}

shared_ptr<Expression> operator/(shared_ptr<Expression> lhs, shared_ptr<Expression> rhs) {
    if (lhs->to_str() == rhs->to_str()) {
        return make_shared<Scalar>(1);
    }
    if (rhs->is_scalar()) {
        auto rhs_scalar = dynamic_pointer_cast<Scalar>(rhs);
        return make_shared<Scalar>(1/rhs_scalar->value) * lhs ;
    } if (rhs->is_atom()) {
        if (rhs->is_reciprocal()) {
            auto rhs_rec = dynamic_pointer_cast<Reciprocal>(rhs);
            return lhs * rhs_rec->content;
        }
        return lhs * make_shared<Reciprocal>(rhs);
    } if (rhs->is_product()) {
        auto temp = make_shared<Multiplication>();
        auto rhs_prod = dynamic_pointer_cast<Multiplication>(rhs);
        temp->scalar = 1 / rhs_prod->scalar;
        for (auto factor : rhs_prod->factors) {
            if (factor->is_reciprocal()) {
                auto factor_rec = dynamic_pointer_cast<Reciprocal>(factor);
                temp->factors.emplace_back(factor_rec->content);
            } else {
                temp->factors.emplace_back(make_shared<Reciprocal>(factor));
            }
        }
        return lhs * temp;
    }
    return lhs * make_shared<Reciprocal>(rhs);
}

shared_ptr<Expression> operator/(double lhs, shared_ptr<Expression> other) {
    return make_shared<Scalar>(lhs) / other;
}

shared_ptr<Expression> operator/(shared_ptr<Expression> lhs, double other) {
    return lhs / make_shared<Scalar>(other);
}


shared_ptr<Expression> operator+(shared_ptr<Expression> lhs, shared_ptr<Expression> rhs) {
    auto ans = make_shared<Sum>();
    if (lhs->is_atom() || lhs->is_product() || lhs->is_scalar()){
        ans->add_sumand(lhs);
    } else {
        auto lhs_sum = dynamic_pointer_cast<Sum>(lhs);
        ans->independent += lhs_sum->independent;
        for (auto sumand : lhs_sum->sumands) {
            ans->add_sumand(sumand);
        }
    }
    if (rhs->is_atom() || rhs->is_product() || rhs->is_scalar()){
        ans->add_sumand(rhs);
    } else {
        auto rhs_sum = dynamic_pointer_cast<Sum>(rhs);
        ans->independent += rhs_sum->independent;
        for (auto sumand : rhs_sum->sumands) {
            ans->add_sumand(sumand);
        }
    }
    if (is_zero(ans->independent) && ans->sumands.size()==1) {
        return ans->sumands[0];
    } if (ans->sumands.size() == 0) {
        return make_shared<Scalar>(ans->independent);
    }
    return ans;
}

shared_ptr<Expression> operator+(double lhs, shared_ptr<Expression> other) {
    return make_shared<Scalar>(lhs) + other;
}

shared_ptr<Expression> operator+(shared_ptr<Expression> lhs, double other) {
    return make_shared<Scalar>(other) + lhs;
}


shared_ptr<Expression> operator-(shared_ptr<Expression> lhs) {
    return make_shared<Scalar>(-1) * lhs;
}

shared_ptr<Expression> operator-(shared_ptr<Expression> lhs, shared_ptr<Expression> rhs) {
    return lhs + (-rhs);
}

shared_ptr<Expression> operator-(double lhs, shared_ptr<Expression> other) {
    return make_shared<Scalar>(lhs) - other;
}

shared_ptr<Expression> operator-(shared_ptr<Expression> lhs, double other) {
    return lhs - make_shared<Scalar>(other);
}

string to_string(shared_ptr<Expression> ex) {
    return ex->to_str();
}

double abs(shared_ptr<Expression> ex) {
    return ex->to_abs();
}

bool is_zero(shared_ptr<Expression> ex) {
    return ex->isZero();
}


