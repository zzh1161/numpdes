#ifndef __GENERIC_FACTORY_HPP__
#define __GENERIC_FACTORY_HPP__

#include <map>
#include <memory>
#include <stdexcept>

template <class AbstractProduct, typename IdType,
            typename ProductCreator = std::unique_ptr<AbstractProduct>(*)()>
class GenericFactory
{
private:
    using CallbackMap = std::map<IdType, ProductCreator>;

public:
    static GenericFactory &createFactory()
    {
        static GenericFactory object;
        return object;
    }
    //return true if registration was successful
    bool registerProduct(IdType Id, ProductCreator createFn)
    {
        return callbacks_.insert(typename CallbackMap::value_type(Id, createFn)).second;
    }
    //return true if the Id was registered before
    bool unregisterProduct(IdType Id)
    {
        return callbacks_.erase(Id) == 1;
    }

    //support construction
    template <class... TS>
    std::unique_ptr<AbstractProduct> createProduct(IdType Id, TS &&...args)
    {
        auto it = callbacks_.find(Id);
        if (it == callbacks_.end())
        {
            throw std::runtime_error("Unknown product ID. ");
        }
        return (it->second)(std::forward<TS>(args)...);
    }

private:
    GenericFactory() = default;
    GenericFactory(const GenericFactory &) = default;
    GenericFactory &operator=(const GenericFactory &) = default;
    ~GenericFactory() = default;

private:
    CallbackMap callbacks_;
};

#endif